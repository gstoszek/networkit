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
#include <algorithm>
#include "EffectiveResistanceDistance.h"
#include "ERDLevel.h"
#include "../components/ConnectedComponents.h"
#include <armadillo>

namespace NetworKit {

   CurrentFlowGroupCloseness::CurrentFlowGroupCloseness(const Graph& G,const count k, const count CB,const double epsilon) : Centrality(G, true),k(k),CB(CB),epsilon(epsilon){
     if (G.isDirected()) throw std::runtime_error("Graph is directed graphs!");
     ConnectedComponents cc(G);
     cc.run();
     if (cc.getPartition().numberOfSubsets() > 1) throw std::runtime_error("Graph has more then one component!");
     if(k>=G.numberOfNodes()) throw std::runtime_error("Size of Group greater then number of nodes!");

     S.resize(k);
     n=G.numberOfNodes();
     TopMatch.resize(n);
     L.set_size(n,n);
     L.zeros();

     vecOfNodes=G.nodes();
     ERD.M=L;
     for(count i=0;i<L.n_rows;i++){
       for(count j=i+1;j<L.n_rows;j++){
         if(G.hasEdge(vecOfNodes[i],vecOfNodes[j])){
           L(i,j)=-1.;
           L(j,i)=-1.;
           L(i,i)+=1;
           L(j,j)+=1;
           TopMatch[i].push_back(j);
           TopMatch[j].push_back(i);
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
      if(CB>0){

          minDegree=1;
          std::vector<std::pair<count,count>> c_indices;
          std::vector<std::pair<node,node>> Matching;
          ID=0;
          c_indices.resize(0);

          ERDLevel Level(ID,vecOfNodes,Matching,TopMatch);
          LevelList.push_back (Level);
          /*update*/
          minDegree=updateMinDegree(minDegree);

          c_indices=peripheralCoarsingIndices();
          Matching=updateMatching(c_indices);
          coarseLaplacian(c_indices);
          ID++;
          Level.set(ID,vecOfNodes,Matching,TopMatch);
          LevelList.push_back (Level);
          minDegree=updateMinDegree(minDegree);
          TopMatch=updateTopMatch(minDegree);

          while((minDegree==1)&&(CB>1)){
            c_indices=peripheralCoarsingIndices();
            Matching=updateMatching(c_indices);
            coarseLaplacian(c_indices);
            ID++;
            Level.set(ID,vecOfNodes,Matching,TopMatch);
            LevelList.push_back (Level);
            minDegree=updateMinDegree(minDegree);
            TopMatch=updateTopMatch(minDegree);
          }
          while(minDegree<CB){
            c_indices=coarsingIndices(minDegree, false);
            Matching=updateMatching(c_indices);
            coarseLaplacian(c_indices);
            ID++;
            Level.set(ID,vecOfNodes,Matching,TopMatch);
            LevelList.push_back (Level);
            minDegree=updateMinDegree(minDegree);
            TopMatch=updateTopMatch(minDegree);
          }
        }
        L=arma::pinv(L,0.01);
        ERD.computeFromPinvL(L,vecOfNodes);
        end = std::chrono::high_resolution_clock::now();
        diff = end-start;
        std::cout << "Initial EffectiveResistanceDistanceMatrice finished in " << diff.count() << "(s)" << "\n\n";
        if(CB>1){
          ID=LevelList[LevelList.size()-1].get_ID();
          while(ID>1){
            Matching=LevelList[ID].get_Matching();
            for(count i=0;i<Matching.size();i++){
              uncoarseEfffectiveResistanceDistanceMatrix(Matching[i].second, Matching[i].first, ID);
            }
            ID--;
          }
          nPeripheralMerges=mergePeripheralNodes();
        }
        end = std::chrono::high_resolution_clock::now();
        diff = end-start;
        std::cout << "Computation of EffectiveResistanceDistanceMatrice finished in " << diff.count() << "(s)" << "\n\n";

        greedy(nPeripheralMerges);
      }
      std::vector<node> CurrentFlowGroupCloseness::getNodesofGroup(){
        return S;
      }
      double CurrentFlowGroupCloseness::getCFGCC() {
        return CFGCC;
      }
      void CurrentFlowGroupCloseness::greedy(count nPeripheralMerges){
        count sampleSize;
        node v;
        node w;
        node s;
        double centrality;
        std::vector<bool> V;
        std::vector<double> mindst;
        std::vector<double> dst;
        std::vector<double> bst;
        std::vector<double> zeroVec;
        std::vector<node> vecOfSamples;

        vecOfSamples = vecOfNodes;
        vecOfSamples.erase(vecOfSamples.end()-nPeripheralMerges, vecOfSamples.end());
        sampleSize=(count)(log(vecOfSamples.size())/(2*epsilon*epsilon));

        if(sampleSize>vecOfSamples.size()){sampleSize=vecOfSamples.size();}

        CFGCC=ERD.M.max()*n;
        V.resize(G.numberOfNodes(),true);
        mindst.resize(G.numberOfNodes(),ERD.M.max());
        zeroVec.resize(G.numberOfNodes(),0.);

        for(count i=0;i<k;i++){
          for (count j=0; j<vecOfSamples.size();j++) {
            if(V[vecOfNodes[j]]){
              std::random_shuffle (vecOfSamples.begin(), vecOfSamples.end());
              v=vecOfNodes[j];
              dst=mindst;
              centrality = 0.;
              for (count l = 0; l < sampleSize; l++) {
                w=vecOfNodes[l];
                if (ERD.M(v,w)< mindst[w])
                  dst[w]=ERD.M(v,w);
                centrality +=dst[w];
              }
              centrality*=((double)(vecOfSamples.size())/(double)(sampleSize));
              for (count l = 0 ; l < nPeripheralMerges; l++) {
                w=vecOfNodes[vecOfNodes.size()-nPeripheralMerges+l];
                if (ERD.M(v,w)< mindst[w])
                  dst[w]=ERD.M(v,w);
                centrality += dst[w];
              }
              if (centrality < CFGCC) {
                CFGCC = centrality;
                bst=dst;
                s = v;
              }
            }
          }
          S[i]=s;
          V[s]=false;
          mindst=bst;
        }
        CFGCC = (double)(n)/CFGCC;
      }

    std::vector<std::vector<node>> CurrentFlowGroupCloseness::updateTopMatch(count minDegree){
      bool search;
      count j;

      node v;
      node w;
      std::vector<std::vector<node>> update_List;
      update_List.resize(G.numberOfNodes());
      for(count i=0;i<vecOfNodes.size();i++){
        v=vecOfNodes[i];
        if(L(i,i)==minDegree){
          search=true;
          j=0;
          while( (j<vecOfNodes.size())&&(search)){
            if((L(i,j)!=0) && (i!=j)){
              w=vecOfNodes[j];
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
      if(i==0){
      std::cout<<"Matrix aufgefressen!" << "\n\n";
      }
      if(min==0){
        std::cout<<"NETWORK_ERROR: min= "<< min << "\n\n";
      }
      return min;
    }
    /***************************************************************************/
    std::vector<std::pair<node,node>> CurrentFlowGroupCloseness::updateMatching(std::vector<std::pair<count,count>> indices){
      std::vector<std::pair<node,node>> Matching;
      Matching.resize(indices.size());
        for(count i=0;i<Matching.size();i++){
          Matching[i].first=vecOfNodes[indices[i].first];
          Matching[i].second=vecOfNodes[indices[i].second];
        }
      return Matching;
    }
    /***************************************************************************/
    std::vector<std::pair<count,count>> CurrentFlowGroupCloseness::peripheralCoarsingIndices(){
      count v;
      count s;

      std::vector<std::pair<count,count>> indices;
      std::vector<count> reverse;
      reverse.resize(G.numberOfNodes());
      for(count i=0;i<vecOfNodes.size();i++){
        reverse[vecOfNodes[i]]=i;
      }
      indices.resize(0);
      for(count i=0;i<L.n_rows;i++){
        if(L(i,i)==1){
          v=vecOfNodes[i];
          s=reverse[TopMatch[v][0]];
          indices.push_back(std::make_pair(i,s));
        }
      }
      return indices;
    }
    /***************************************************************************/
    std::vector<std::pair<count,count>> CurrentFlowGroupCloseness::coarsingIndices(count cDegree, bool Random){
        bool search;
        count c_index;
        count s_index;
        count l;
        node c;
        node s;
        std::vector<bool> s_List;
        std::vector<count> c_List;
        std::vector<std::pair<count,count>> indices;
        std::vector<count> reverse;

        c_List.resize(0);
        s_List.resize(G.numberOfNodes(),true);
        reverse.resize(G.numberOfNodes());
        for(count i=0;i<vecOfNodes.size();i++){
          reverse[vecOfNodes[i]]=i;
        }
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
          c=vecOfNodes[c_index];
          if(s_List[c]){
            search=true;
            l=0;
            while((search)&&(l<TopMatch[c].size())){
              s=TopMatch[c][l];
              s_index=reverse[s];
              if(s_List[s]){
                s_List[s]=false;
                search=false;
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
    void CurrentFlowGroupCloseness::uncoarseEfffectiveResistanceDistanceMatrix(node u,node v,count ID){
      bool search;
      count l;
      node w;
      std::vector<node> uAdj;
      std::vector<node> vAdj;

      uAdj= LevelList[ID].getVecOfAdj(u);
      vAdj= LevelList[ID].getVecOfAdj(v);

      ERD.firstJoin(vecOfNodes,u,v,1.);
      vecOfNodes.push_back (v);

      for(count i=1;i<vAdj.size();i++){
        w=vAdj[i];
        search=true;
        l=0;
        while((l<uAdj.size())&&(search)){
         if(uAdj[l]=w){
          ERD.edgeFire(vecOfNodes,v,w,1.);
          search=false;
          }
          else{
            l++;
          }
        }
        if(search==true){
          ERD.corollary(vecOfNodes,u,v,w,1.,1.);
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
      merge_value.resize(G.numberOfNodes(),0);
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
        ERD.firstJoin(vecOfNodes,s,v,merge_value[s]);
        vecOfNodes.push_back (v);
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
      reverse.resize(G.numberOfNodes());

      for(count i=0;i<vecOfNodes.size();i++){
        reverse[vecOfNodes[i]]=i;
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
        }
      }
      arma::uvec indices(vecOfNodes.size()-Matchings.size());
      a=0;
      l=0;
      for(count i=0;i<vecOfNodes.size();i++){
        if(Matchings[l].first==i){
          if(l<Matchings.size()-1)
            l++;
        }
        else{
          indices(a)=i;
          a++;
        }
      }
      for(count i=0;i<Matchings.size();i++){
        vecOfNodes.erase(vecOfNodes.begin()+Matchings[i].first-i);
      }
      L=L.submat(indices, indices);

    }
} /* namespace NetworKit*/
