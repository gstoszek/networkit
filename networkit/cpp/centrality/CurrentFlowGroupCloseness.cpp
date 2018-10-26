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
#include "ERDLevel.h"
#include <armadillo>

namespace NetworKit {

   CurrentFlowGroupCloseness::CurrentFlowGroupCloseness(const Graph& G,const count k, const count CB) : Centrality(G, true),k(k),CB(CB){
        S.resize(k);
        CFGCC = 0.;
        n=G.upperNodeIdBound();
        vList.resize(n);
        TopMatch.resize(n);
        /*Laplacian*/
        L.set_size(n,n);
        L.zeros();

        ERD.M=L;
        Adj=L;

        for(count i=0;i<L.n_rows;i++){
          vList[i]=i;
          for(count j=0;j<L.n_rows;j++){
            if(G.hasEdge(i,j)){
              L(i,j)=-1.;
              L(i,i)+=1;
              TopMatch[i].push_back(j);
            }
          }
        }
    }

    void CurrentFlowGroupCloseness::run() {
      count ID;
      count nPeripheralMerges;
      count minDegree;
      auto start = std::chrono::high_resolution_clock::now();
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff;
      std::vector<std::pair<node,node>> Matching;
      std::vector<std::pair<count,count>> Indices;

      nPeripheralMerges=0;
      if(CB>1){
          minDegree=1;
          std::vector<std::pair<count,count>> c_indices;
          std::vector<std::pair<node,node>> Matching;
          ID=0;
          c_indices.resize(0);

          ERDLevel Level(ID,vList,Matching,TopMatch,Adj);
          LevelList.push_back (Level);
          /*update*/
          minDegree=updateMinDegree(minDegree);

          while(minDegree==1){
            c_indices=peripheralCoarsingIndices();
            Matching=updateMatching(c_indices);
            coarseLaplacian(c_indices);
            ID++;
            Level.set(ID,vList,Matching,TopMatch,Adj);
            LevelList.push_back (Level);
            minDegree=updateMinDegree(minDegree);
            TopMatch=updateTopMatch(minDegree);
          }
          /*
          while(minDegree<CB){
            c_indices=coarsingIndices(minDegree, false);
            Matching=updateMatching(c_indices);
            coarseLaplacian(c_indices);
            ID++;
            Level.set(ID,vList,Matching,TopMatch,Adj);
            LevelList.push_back (Level);
            minDegree=updateMinDegree(minDegree);
            TopMatch=updateTopMatch(minDegree);
          }
          */
        }
        L=arma::pinv(L, 0.01);
        ERD.computeFromPinvL(L,vList);
        if(CB>1){
          ID=LevelList[LevelList.size()-1].get_ID();
          while(ID>1){
            Matching=LevelList[ID].get_Matching();
            for(count i=0;i<Matching.size();i++){
              uncoarse(Matching[i].second, Matching[i].first, ID);
            }
            ID--;
            vList=LevelList[ID].get_vList();
          }
          nPeripheralMerges=mergePeripheralNodes();
        }
        end = std::chrono::high_resolution_clock::now();
        diff = end-start;
        std::cout << "Computation of EffectiveResistanceDistanceMatrice finished in " << diff.count() << "(s)" << "\n\n";


        greedy(nPeripheralMerges);
      }
      /*******************************************************************************************************************************************************/
      std::vector<node> CurrentFlowGroupCloseness::getNodesofGroup(){
          return S;
      }
      /*******************************************************************************************************************************************************/
      double CurrentFlowGroupCloseness::getCFGCC() {
          return CFGCC;
      }
      /**************************************************************************************************************/
        void CurrentFlowGroupCloseness::cleanNetwork(){
          count n2;
          n2=0;
          for(count i=0;i<n;i++){
            if(L(i,i)==0){
              vList.erase(vList.begin()+i-n2);
              n2++;
            }
          }
          if(n2>0){
            n-=n2;
            arma::uvec indices(n);
            for(count i=0;i<n;i++){
                indices(i)=vList[i];
            }
            L=L.submat(indices, indices);
          }
        }
        /**************************************************************************************************************/
        void CurrentFlowGroupCloseness::greedy(count n_peripheral_merges){
          count k_max;
          node s;
          node s_next;
          node v;
          double S_CFGCC;
          double S_currentCFGCC;
          double scaling_factor;
          std::vector<bool> V;
          std::vector<double> d;

          k_max=S.size();
          scaling_factor=(double)(n);
          if(k_max>vList.size()){
            k_max=vList.size()-1;
            scaling_factor=(double)((G.upperNodeIdBound()))/(double)(k_max);
          }
          S_CFGCC=0.;
          V.resize(vList.size(),true);
          d.resize(G.upperNodeIdBound(),n*n);
          for(count i=0;i<k_max;i++){
            /*Maximal Gain Loop*/
            for (count j=0; j<vList.size()-n_peripheral_merges;j++) {
                //use vector of bools
                if (V[j]){
                    s=vList[j];
                    /*Sample Loop*/
                    S_currentCFGCC = 0.;
                    for (count l = 0; l < vList.size(); l++) {
                        v=vList[l];
                        if (ERD.M(v,s)< d[v])
                            S_currentCFGCC = S_currentCFGCC + ERD.M(v,s);
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
                if (ERD.M(v,s_next)< d[v]) {
                    d[v] = ERD.M(v,s_next);
                }
            }
            S[i]=s_next;
            V[s_next]=false;
        }
        CFGCC=S_CFGCC;
    }

    void CurrentFlowGroupCloseness::computeInitialERD(count CB){
      count ID;
      count minDegree;

      minDegree=1;
      std::vector<std::pair<count,count>> c_indices;
      std::vector<std::pair<node,node>> Matching;
      ID=0;
      c_indices.resize(0);

      ERDLevel Level(ID,vList,Matching,TopMatch,Adj);
      LevelList.push_back (Level);
      /*update*/
      minDegree=updateMinDegree(minDegree);

      while(minDegree==1){
        c_indices=peripheralCoarsingIndices();
        Matching=updateMatching(c_indices);
        coarseLaplacian(c_indices);
        ID++;
        Level.set(ID,vList,Matching,TopMatch,Adj);
        LevelList.push_back (Level);
        minDegree=updateMinDegree(minDegree);
        TopMatch=updateTopMatch(minDegree);
      }
      /*
      while(minDegree<CB){
        c_indices=coarsingIndices(minDegree, false);
        Matching=updateMatching(c_indices);
        coarseLaplacian(c_indices);
        ID++;
        Level.set(ID,vList,Matching,TopMatch,Adj);
        LevelList.push_back (Level);
        minDegree=updateMinDegree(minDegree);
        TopMatch=updateTopMatch(minDegree);
      }
      */
    }
    /***************************************************************************/
    std::vector<std::vector<node>> CurrentFlowGroupCloseness::updateTopMatch(count minDegree){
      bool search;
      count j;

      node v;
      node w;
      std::vector<std::vector<node>> update_List;
      update_List.resize(G.upperNodeIdBound());
      for(count i=0;i<vList.size();i++){
        v=vList[i];
        if(L(i,i)==minDegree){
          search=true;
          j=0;
          while( (j<vList.size())&&(search)){
            if((L(i,j)!=0) && (i!=j)){
              w=vList[j];
              update_List[v].push_back(w);
              if(update_List[v].size()==minDegree){
                search=false;
              }
            }
            j++;
          }
        }
      }
      return update_List;
    }
    /***************************************************************************/
    count CurrentFlowGroupCloseness::updateMinDegree(count minDegree){
      bool search;
      count min;
      count i;

      search=true;
      min=n;
      i=0;

      while((i<L.n_rows)&&(search)){
        if(L(i,i)<min){
          min=L(i,i);
          if(min==minDegree){
            search=false;
          }
        }
        i++;
      }
      if(min==0){
        std::cout<<"NETWORK_ERROR: min= "<< min << "\n\n";
      }
      for(count i=0;i<L.n_rows;i++){
        if(L(i,i)==0){
          std::cout<<"ERROR\n";
        }
      }
      return min;
    }
    /***************************************************************************/
    std::vector<std::pair<node,node>> CurrentFlowGroupCloseness::updateMatching(std::vector<std::pair<count,count>> indices){
      std::vector<std::pair<node,node>> Matching;
      Matching.resize(indices.size());
        for(count i=0;i<Matching.size();i++){
          Matching[i].first=vList[indices[i].first];
          Matching[i].second=vList[indices[i].second];
        }
      return Matching;
    }
    /***************************************************************************/
    std::vector<std::pair<count,count>> CurrentFlowGroupCloseness::peripheralCoarsingIndices(){
      count v;
      count s;

      std::vector<std::pair<count,count>> indices;
      std::vector<count> reverse;
      reverse.resize(G.upperNodeIdBound());
      for(count i=0;i<vList.size();i++){
        reverse[vList[i]]=i;
      }
      indices.resize(0);
      for(count i=0;i<L.n_rows;i++){
        if(L(i,i)==1){
          v=vList[i];
          s=reverse[TopMatch[v][0]];
          indices.push_back(std::make_pair(i,s));
        }
      }
      return indices;
    }
    /***************************************************************************/
    std::vector<std::pair<count,count>> CurrentFlowGroupCloseness::coarsingIndices(count cDegree, bool Random){
        bool s_found;

        count c_index;
        count s_index;
        count l;

        node c;
        /*free supernodes*/
        std::vector<bool> s_List;
        /*potential candidates*/
        std::vector<count> c_List;
        /*c_index-s_index mapping*/
        std::vector<count> reverse;
        reverse.resize(G.upperNodeIdBound());
        for(count i=0;i<vList.size();i++){
          reverse[vList[i]]=i;
        }
        std::vector<std::pair<count,count>> indices;

        c_List.resize(0);
        s_List.resize(G.upperNodeIdBound(),true);
        /*potential candidates*/
        for(count i=0;i<L.n_rows;i++){
          if(L(i,i)==cDegree){
            c_List.push_back(i);
          }
        }
        if(Random){
          std::random_shuffle (c_List.begin(), c_List.end());
        }
        for(count i=0;i<c_List.size();i++){
          c_index=c_List[i];
          c=vList[c_index];
          if(s_List[c_index]){
            s_found=false;
            l=0;
            while(!(s_found)&&(k<TopMatch[c].size())){
              s_index=reverse[TopMatch[c][l]];
              if(s_List[s_index]){
                s_List[s_index]=false;
                s_found=true;
                indices.push_back(std::make_pair(c_index,s_index));
              }
              else{
                l++;
              }
            }
          }
        }
        return indices;
      }
    /***************************************************************************/
    void CurrentFlowGroupCloseness::uncoarse(node s,node v,count ID){
      node w;
      std::vector<node> v_TopMatch;
      arma::vec v_Adj;
      arma::vec s_Adj;

      v_TopMatch = LevelList[ID].get_vecofTopmatch(v);
      v_Adj = LevelList[ID-1].get_CalofAdj(v);
      s_Adj = LevelList[ID-1].get_CalofAdj(s);
      ERD.firstJoin(vList,s,v,1.);
      vList.push_back (v);
      LevelList[ID].set_i_j_ofAdj(s,v,s_Adj(v));

      for(count i=1;i<v_TopMatch.size();i++){
        w=v_TopMatch[i];
        if(v_Adj(w)!=0){
          if(s_Adj(w)!=0){
            ERD.edgeFire(vList,v,w,1.);
            LevelList[ID].set_i_j_ofAdj(v,w,v_Adj(w));
            LevelList[ID].set_i_j_ofAdj(v,w,s_Adj(w)-v_Adj(w));
          }
          else{
            ERD.edgeFire(vList,v,w,1.);
            LevelList[ID].set_i_j_ofAdj(v,w,v_Adj(w));
            ERD.nonBridgeDelete(vList,s,w,1.);
            LevelList[ID].set_i_j_ofAdj(v,w,0.);
          }
        }
      }
    }
    /***************************************************************************/
    count CurrentFlowGroupCloseness::mergePeripheralNodes(){
      node s;
      node v;
      std::vector<std::pair<node,node>> merge;
      std::vector<std::pair<node,node>> Matching;
      std::vector<count> merge_value;

      Matching=LevelList[1].get_Matching();
      merge_value.resize(G.upperNodeIdBound(),0);
      for(count i=0;i<Matching.size();i++){
        v=Matching[i].first;
        s=Matching[i].second;
        merge_value[s]++;
        if(merge_value[s]==1){
          merge.push_back (std::make_pair(v,s));
        }
      }
      for(count i=0;i<merge.size();i++){
        v=merge[i].first;
        s=merge[i].second;
        ERD.firstJoin(vList,s,v,merge_value[s]);
        LevelList[1].set_i_j_ofAdj(s,v,merge_value[s]);
        vList.push_back (v);
      }
      return merge.size();
    }
    /***************************************************************************/
    void CurrentFlowGroupCloseness::coarseLaplacian(std::vector<std::pair<count,count>> Matchings){
      count a;
      count l;
      /*s=supernode - c = to be coarsed node - w = adjacent node to c*/
      node s;
      node c;
      node w;

      std::vector<count> reverse;
      reverse.resize(G.upperNodeIdBound());
      for(count i=0;i<vList.size();i++){
        reverse[vList[i]]=i;
      }
      for(count i=0;i<Matchings.size();i++){
        c=Matchings[i].first;
        s=Matchings[i].second;
        L(s,s)-=1;
        for(count j=1;j<TopMatch[c].size();j++){
          w=reverse[TopMatch[c][j]];
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
      arma::uvec indices(vList.size()-Matchings.size());
      a=0;
      l=0;
      for(count i=0;i<vList.size();i++){
        if(Matchings[l].first==i){
          l++;
        }
        else{
          indices(a)=i;
          a++;
        }
      }
      for(count i=0;i<Matchings.size();i++){
        vList.erase(vList.begin()+Matchings[i].first-i);
      }
      L=L.submat(indices, indices);

    }
} /* namespace NetworKit*/
