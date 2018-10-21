/*
 * CurrentFlowGroupCloseness.cpp
 *
 *      Author: gstoszek
 */

#include "CurrentFlowGroupCloseness.h"
#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include "../auxiliary/Log.h"
#include <chrono>
#include <stdlib.h>
#include <cmath>
#include "EffectiveResistanceDistance.h"
#include "ERD2.h"
#include "ERDLevel.h"
#include <armadillo>

namespace NetworKit {
/*
 * ***SAMPLING*** (Number *** II ***)
 */

   CurrentFlowGroupCloseness::CurrentFlowGroupCloseness(const Graph& G, const double epsilon, const double delta, const count groupsize, const count upperDegreeBound) : Centrality(G, true), epsilon(epsilon), delta(delta), groupsize(groupsize), upperDegreeBound(upperDegreeBound){
        count n;

        S.clear();
        S.resize(groupsize);

        CFGCC = 0.;

        n=G.upperNodeIdBound();
        vList.resize(n);

        /*Laplacian*/
        L.set_size(n,n);
        L.zeros();
        /*EffectiveResistanceDistance Matrix ERD*/
        ERD=L;
        Adj=L;
        vList.resize(n);

        G.forNodes([&](node v){
            vList[v]=v;
            G.forNodes([&](node w){
              if(G.hasEdge(v,w)){
                L(v,w)=-1.;
                L(v,v)+=1.;
                Adj(v,w)=1;
              }
            });
        });
    }
    /*******************************************************************************************************************************************************/
    std::vector<node> CurrentFlowGroupCloseness::getNodesofGroup(){
        return S;
    }
    /*******************************************************************************************************************************************************/
    double CurrentFlowGroupCloseness::getCFGCC() {
        return CFGCC;
    }
    /*******************************************************************************************************************************************************/
    void CurrentFlowGroupCloseness::setCFGCC(double newCFGCC) {
        CFGCC = newCFGCC;
    }
    /*******************************************************************************************************************************************************/
    void CurrentFlowGroupCloseness::run() {

        count n;
        count ID;
        auto start = std::chrono::high_resolution_clock::now();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff;

        std::vector<std::pair<node,node>> Matching;
        arma::vec s_Adj;
        arma::vec v_Adj;
        n = G.upperNodeIdBound();

        std::cout << "Number of Nodes:" << n << "\n";
        /**************************************************************************************************************/
        /*calculating D-Matrix with ERD2*/
        /**************************************************************************************************************/
        std::cout << "Start: Computing initial ERD with the upper Degree Bound: "<<upperDegreeBound<< "\n";
        if(upperDegreeBound<2){
          compute_ERD();
          end = std::chrono::high_resolution_clock::now();
          diff = end-start;
          std::cout << "Finished in " << diff.count() << "(s)" << "\n\n";
          update_CFGCC();
        }
        else{
        compute_initial_ERD(upperDegreeBound);
        end = std::chrono::high_resolution_clock::now();
        diff = end-start;
        std::cout << "Finished in " << diff.count() << "(s)" << "\n\n";

        ID=LevelList[LevelList.size()-1].get_ID();
        while(ID>0){
          Matching=LevelList[ID].get_Matching();
          for(count i=0;i<Matching.size();i++){
            uncoarse(Matching[i].second, Matching[i].first);
          }
          ID--;
          vList=LevelList[ID].get_vList();
          //std::cout << "CFGCC: "<< CFGCC << ", Number of Nodes: " << vList.size()<< "\n";
        }
        update_CFGCC();
        }
      }
        /**************************************************************************************************************/
        /*Greedy*/
        /**************************************************************************************************************/
        void CurrentFlowGroupCloseness::update_CFGCC(){
          count n;
          count k_max;
          node s;
          node s_next;
          node v;
          double S_CFGCC;
          double S_currentCFGCC;
          double scaling_factor;
          std::vector<bool> V;
          std::vector<double> d;

          n=G.upperNodeIdBound();
          k_max=S.size();
          scaling_factor=(double)(n);
          if(k_max>vList.size()){
            k_max=vList.size()-1;
            scaling_factor=(double)(n)/(double)(k_max);
          }
          S_CFGCC=0.;
          V.resize(vList.size(),true);
          d.resize(n,(double) (n*n));
          for(count i=0;i<k_max;i++){
            /*Maximal Gain Loop*/
            for (count j=0; j<vList.size();j++) {
                //use vector of bools
                if (V[j]){
                    s=vList[j];
                    /*Sample Loop*/
                    S_currentCFGCC = 0.;
                    for (count l = 0; l < vList.size(); l++) {
                        v=vList[l];
                        if (ERD(v,s)< d[v])
                            S_currentCFGCC = S_currentCFGCC + ERD(v,s);
                        else {
                            S_currentCFGCC = S_currentCFGCC + d[v];
                        }
                    }
                    S_currentCFGCC =  scaling_factor / S_currentCFGCC;
                    if (S_currentCFGCC > S_CFGCC) {
                        S_CFGCC = S_currentCFGCC;
                        s_next = s;
                    }
                }
            }
            for(count j= 0; j<vList.size(); j++) {
              v=vList[j];
                if (ERD(v,s_next)< d[v]) {
                    d[v] = ERD(v,s_next);
                }
            }
            S[i]=s_next+1;
            V[s_next]=false;
        }
        CFGCC=S_CFGCC;
    }

    void CurrentFlowGroupCloseness::compute_initial_ERD(count upperDegreeBound){
      count ID;
      count minDegree;

      node v;
      node w;

      std::vector<std::pair<count,count>> c_indices;
      std::vector<std::pair<node,node>> Matching;
      ID=0;
      c_indices.resize(0);

     ERDLevel Level(ID,vList,c_indices,Adj);
     LevelList.push_back (Level);
      /*update*/
      minDegree=update_minDegree();
      TopMatch=update_TopMatch();
      while(minDegree==1){
        c_indices=peripheral_indices();
        Matching=update_Matching(c_indices);
        coarse_L(c_indices);
        ID++;
        Level.set_ID(ID);
        Level.set_vList(vList);
        Level.set_Matching(Matching);
        Level.set_Adj(Adj);
        LevelList.push_back (Level);
        minDegree=update_minDegree();
        TopMatch=update_TopMatch();
        std::cout <<"Round:"<< ID <<"\n";
      }
      while(minDegree<upperDegreeBound){
        c_indices=coarsing_indices(minDegree, true);
        Matching=update_Matching(c_indices);
        coarse_L(c_indices);
        /*update*/
        ID++;
        Level.set_ID(ID);
        Level.set_vList(vList);
        Level.set_Matching(Matching);
        Level.set_Adj(Adj);
        LevelList.push_back (Level);
        TopMatch=update_TopMatch();
        minDegree=update_minDegree();
      }
      /*create pseudo inverse*/
      L=arma::pinv(L, 0.01);
      /*calculate initial ERD Matrix*/
      for(count i=0;i<vList.size();i++){
        v=vList[i];
        for(count j=i+1;j<vList.size();j++){
          w=vList[j];
          ERD(v,w)=L(i,i)+L(j,j)-2.*L(i,j);
          ERD(w,v)=ERD(v,w);
        }
      }
    }
    /***************************************************************************/
    void CurrentFlowGroupCloseness::compute_ERD(){
      node v;
      node w;
      L=arma::pinv(L, 0.01);
      /*calculate initial ERD Matrix*/
      for(count i=0;i<vList.size();i++){
        v=vList[i];
        for(count j=i+1;j<vList.size();j++){
          w=vList[j];
          ERD(v,w)=L(i,i)+L(j,j)-2.*L(i,j);
          ERD(w,v)=ERD(v,w);
        }
      }
    }
    /***************************************************************************/
    std::vector<std::vector<node>> CurrentFlowGroupCloseness::update_TopMatch(){
      std::vector<std::vector<node>> update_List;
      update_List.resize(vList.size());
      for(count i=0;i<L.n_rows;i++){
        for(count j=i+1;j<L.n_rows;j++){
          if(L(i,j)!=0){
            update_List[j].push_back(i);
            update_List[i].push_back(j);
          }
        }
      }
      return update_List;
    }
    /***************************************************************************/
    count CurrentFlowGroupCloseness::update_minDegree(){
      count min;
      min=G.upperNodeIdBound();
      for(count i=0;i<L.n_rows;i++){
        if(L(i,i)<min){
          min=L(i,i);
        }
      }
      if(min==0){
        std::cout<<"NETWORK_ERROR: min= "<< min << "\n\n";
      }
      return min;
    }
    /***************************************************************************/
    std::vector<std::pair<node,node>> CurrentFlowGroupCloseness::update_Matching(std::vector<std::pair<count,count>> indices){
      std::vector<std::pair<node,node>> Matching;
      Matching.resize(indices.size());
        for(count i=0;i<Matching.size();i++){
          Matching[i].first=vList[indices[i].first];
          Matching[i].second=vList[indices[i].second];
        }
      return Matching;
    }
    /***************************************************************************/
    std::vector<std::pair<count,count>> CurrentFlowGroupCloseness::peripheral_indices(){
      std::vector<std::pair<count,count>> indices;
      indices.resize(0);
      for(count i=0;i<L.n_rows;i++){
        if(L(i,i)==1){
          indices.push_back(std::make_pair(i,TopMatch[i][0]));
        }
      }
      return indices;
    }
    /***************************************************************************/
      /*
      Type 0: Matching without Permutation -> not randomized;
      Type 1: Matching with Permutation -> randomized;
      */
    std::vector<std::pair<count,count>> CurrentFlowGroupCloseness::coarsing_indices(count cDegree, bool Random){
        bool s_found;

        count c_index;
        count s_index;
        count k;
        count ID;

        /*free supernodes*/
        std::vector<bool> s_List;
        /*potential candidates*/
        std::vector<count> c_List;
        /*c_index-s_index mapping*/
        std::vector<std::pair<count,count>> indices;

        ID=0;
        ERDLevel Level(ID,vList,indices,Adj);

        c_List.resize(0);
        s_List.resize(vList.size(),true);
        /*potential candidates*/
        for(count i=0;i<L.n_rows;i++){
          if(L(i,i)==cDegree){
            c_List.push_back(i);
          }
        }
        if(Random){
          std::random_shuffle (c_List.begin(), c_List.end());
        }
        for(count i=0;i<L.n_rows;i++){
          c_index=c_List[i];
          if(s_List[c_index]){
            s_found=false;
            k=0;
            while(!(s_found)&&(k<TopMatch[c_index].size())){
              s_index=TopMatch[c_index][k];
              if(s_List[s_index]){
                s_List[s_index]=false;
                s_found=true;
                indices.push_back(std::make_pair(c_index,s_index));
              }
            }
          }
        }
        return indices;
      }
    /***************************************************************************/
    void CurrentFlowGroupCloseness::uncoarse(node s,node v){
        first_join(s,v);
        vList.push_back (v);
    }
    /***************************************************************************/
    void CurrentFlowGroupCloseness::coarse_L(std::vector<std::pair<count,count>> Matchings){
      count k;
      count l;

      /*s=supernode - c = to be coarsed node - w = adjacent node to c*/
      node s;
      node c;
      node w;

      arma::uvec indices(vList.size()- Matchings.size());
      for(count i=0;i<Matchings.size();i++){
        c=Matchings[i].first;
        s=Matchings[i].second;
        L(s,c)=0;
        L(c,s)=0;
        for(count j=1;j<TopMatch[c].size();j++){
          w=TopMatch[c][j];
          if(L(s,w)==0){
            L(s,w)=-1.;
            L(w,s)=-1.;
            /*Diagonal*/
            L(s,s)+=1.;
          }
          Adj(vList[s],vList[w])+=Adj(vList[c],vList[w]);
          Adj(vList[w],vList[s])=Adj(vList[s],vList[w]);
        }
        Adj(vList[s],vList[c])=0;
        Adj(vList[c],vList[s])=0;
      }
      k=0;
      l=0;
      for(count i=0;i<vList.size();i++){
        if(Matchings[l].first==i){
          l++;
        }
        else{
          indices(k)=i;
          k++;
        }
      }
      for(count i=0;i<Matchings.size();i++){
        vList.erase(vList.begin()+Matchings[i].first-i);
      }
      L=L.submat(indices, indices);
    }
    /**************************************************************************/
    void CurrentFlowGroupCloseness::first_join(node s, node v){
      node w;
      for(count i=0;i<vList.size();i++){
        w=vList[i];
        ERD(w,v)=ERD(w,s)+1.;
        ERD(v,w)=ERD(w,v);
      }
    }
    /**************************************************************************/
    void CurrentFlowGroupCloseness::edge_fire(node v, node w){
      node x;
      node y;
      double fix_factor;
      arma::Mat<double> ERD2;

      fix_factor=4.*(1.+ERD(v,w));
      ERD2 = ERD;

      for(count i=0;i<vList.size();i++){
        x=vList[i];
        for(count j=i+1;j<vList.size();j++){
          y=vList[j];
          ERD2(x,y)=ERD(x,w)-ERD(x,v);
          ERD2(x,y)-=ERD(w,y)-ERD(v,y);
          ERD2(x,y)*=ERD2(x,y);
          ERD2(x,y)/=fix_factor;
          ERD2(x,y)=ERD(x,y)-ERD2(x,y);
          ERD2(y,x)=ERD2(x,y);
        }
      }
      ERD=ERD2;
    }
} /* namespace NetworKit*/
