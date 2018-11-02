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
     Adj.resize(n);
     L.set_size(n,n);
     L.zeros();
     std::vector<node> vecOfNodesMapping;
     vecOfNodes.resize(n);
     vecOfNodesMapping=G.nodes();
     ERD.M=L;
     for(count i=0;i<L.n_rows;i++){
       vecOfNodes[i]=i;
       for(count j=i+1;j<L.n_rows;j++){
         if(G.hasEdge(vecOfNodesMapping[i],vecOfNodesMapping[j])){
           L(i,j)=-1.;
           L(j,i)=-1.;
           L(i,i)+=1;
           L(j,j)+=1;
           Adj[i].push_back(j);
           Adj[j].push_back(i);
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
      std::vector<std::pair<count,count>> c_indices;;
      nPeripheralMerges=0;
      ID=0;
      if(CB>0){
          minDegree=1;
          ID++;
          c_indices=peripheralCoarsingIndices();
          coarseLaplacian(c_indices,ID);
          minDegree=updateMinDegree(1);
          if(CB>1){
            while((minDegree<CB)&&(vecOfNodes.size()>k)){
              ID++;
              if(minDegree==1){
                c_indices=peripheralCoarsingIndices();
              }
              else{
                c_indices=coarsingIndices(minDegree, false);
              }
              coarseLaplacian(c_indices,ID);
              minDegree=updateMinDegree(1);
          }
        }
      }
      L=arma::pinv(L,0.01);
      ERD.computeFromPinvL(L,vecOfNodes);
      end = std::chrono::high_resolution_clock::now();
      diff = end-start;
      std::cout << "Initial EffectiveResistanceDistanceMatrix finished in " << diff.count() << "(s)" << "\n\n";
      if(CB>1){
        ID=LevelList[LevelList.size()-1].getID();
        while(ID>1){
          uncoarseEfffectiveResistanceDistanceMatrix(ID);
          ID--;
        }
        nPeripheralMerges=mergePeripheralNodes();
      }
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
              if (ERD.M(v,w)< mindst[w]){
                dst[w]=ERD.M(v,w);
              }
              centrality +=dst[w];
            }
            centrality*=((double)(vecOfSamples.size())/(double)(sampleSize));
            for (count l = 0 ; l < nPeripheralMerges; l++) {
              w=vecOfNodes[vecOfNodes.size()-nPeripheralMerges+l];
              if (ERD.M(v,w)< mindst[w]){
                dst[w]=ERD.M(v,w);
              }
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
      return min;
    }
    std::vector<std::pair<count,count>> CurrentFlowGroupCloseness::peripheralCoarsingIndices(){
      count s_i;
      node c;

      std::vector<std::pair<count,count>> indices;
      std::vector<count> reverse;
      reverse.resize(G.numberOfNodes());
      for(count i=0;i<vecOfNodes.size();i++){
        reverse[vecOfNodes[i]]=i;
      }
      indices.resize(0);
      for(count c_i=0;c_i<L.n_rows;c_i++){
        if(L(c_i,c_i)==1){
          c=vecOfNodes[c_i];
          s_i=reverse[Adj[c][0]];
          indices.push_back(std::make_pair(c_i,s_i));
        }
      }
      return indices;
    }
    /***************************************************************************/
    std::vector<std::pair<count,count>> CurrentFlowGroupCloseness::coarsingIndices(count cDegree, bool Random){
        bool search;
        count c_i;
        count s_i;
        count l;
        node c;
        node s;
        std::vector<bool> vecOfSupernodes;
        std::vector<count> vecOfCoarseNodes;
        std::vector<std::pair<count,count>> indices;
        std::vector<count> reverse;

        vecOfCoarseNodes.resize(0);
        vecOfSupernodes.resize(G.numberOfNodes(),true);
        reverse.resize(G.numberOfNodes());

        for(count i=0;i<vecOfNodes.size();i++){
          reverse[vecOfNodes[i]]=i;
        }
        for(count i=0;i<L.n_rows;i++){
          if(L(i,i)==cDegree){
            vecOfCoarseNodes.push_back(i);
          }
        }
        if(Random){
          std::random_shuffle (vecOfCoarseNodes.begin(), vecOfCoarseNodes.end());
        }
        for(count i=0;i<vecOfCoarseNodes.size();i++){
          c_i=vecOfCoarseNodes[i];
          c=vecOfNodes[c_i];
          if(vecOfSupernodes[c]){
            search=true;
            l=0;
            while((search)&&(l<Adj[c].size())){
              s=Adj[c][l];
              s_i=reverse[s];
              if(vecOfSupernodes[s]){
                vecOfSupernodes[s]=false;
                vecOfSupernodes[c]=false;
                search=false;
                indices.push_back(std::make_pair(c_i,s_i));
              }
              else{
                l++;
              }
            }//End while
          }
        }
        return indices;
      }
    /***************************************************************************/
    void CurrentFlowGroupCloseness::uncoarseEfffectiveResistanceDistanceMatrix(count ID){
      bool search;
      count l;
      count j;
      node u;
      node v;
      node w;
      std::vector<std::tuple<node,node,count>> vecOfTransferredEdges;
      vecOfTransferredEdges=LevelList[ID-1].getVecOfTransferredEdges();
      for(count i=0;i<vecOfTransferredEdges.size();i++){
        j=vecOfTransferredEdges.size()-1-i;
        if(std::get<2>(vecOfTransferredEdges[j]) == 1){
          u=std::get<1>(vecOfTransferredEdges[j]);
          v=std::get<0>(vecOfTransferredEdges[j]);
          ERD.firstJoin(vecOfNodes,u,v,1.);
          vecOfNodes.push_back (v);
          Adj[v].push_back (u);
        }
        else if(std::get<2>(vecOfTransferredEdges[j])== 2){
          w=std::get<1>(vecOfTransferredEdges[j]);
          ERD.edgeFire(vecOfNodes,v,w,1.);
          Adj[v].push_back (w);
        }
        else if(std::get<2>(vecOfTransferredEdges[j])== 3){
          w=std::get<1>(vecOfTransferredEdges[j]);
          ERD.edgeFire(vecOfNodes,v,w,1.);
          ERD.nonBridgeDelete(vecOfNodes,u,w,1.);
          //ERD.corollary(vecOfNodes,u,v,w,1.,1.);
          Adj[v].push_back (w);
          search = true;
          l= 0;
          while((l<Adj[u].size())&&(search)){
            if(Adj[u][l]==w){
              Adj[u].erase(Adj[u].begin()+l);
              search=false;
            }
            else{
              l++;
            }
          }//end while
        }//end case 3
      }
    }
    /***************************************************************************/
    count CurrentFlowGroupCloseness::mergePeripheralNodes(){
      node s;
      node v;
      std::vector<std::pair<node,node>> merge;
      std::vector<std::tuple<node,node,count>> matchings;
      std::vector<count> merge_value;

      matchings=LevelList[0].getVecOfTransferredEdges();
      merge_value.resize(G.numberOfNodes(),0);
      for(count i=0;i<matchings.size();i++){
        v=std::get<0>(matchings[i]);
        s=std::get<1>(matchings[i]);
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
    void CurrentFlowGroupCloseness::coarseLaplacian(std::vector<std::pair<count,count>> matchings,count ID){
      bool search;
      count a;
      count l;
      count s_i;
      count c_i;
      count w_i;
      /*s=supernode - c = to be coarsed node - w = adjacent node to c*/
      node s;
      node c;
      node w;
      std::vector<bool> reverseBool;
      std::vector<count> reverse;
      std::vector<std::tuple<node,node,count>> vecOfTransferredEdges;

      reverseBool.resize(G.numberOfNodes(),false);
      reverse.resize(G.numberOfNodes());
      vecOfTransferredEdges.resize(0);
      for(count i=0;i<vecOfNodes.size();i++){
        reverse[vecOfNodes[i]]=i;
        reverseBool[vecOfNodes[i]]=true;
      }
      for(count i=0;i<matchings.size();i++){
        c_i=matchings[i].first;
        c=vecOfNodes[c_i];
        s_i=matchings[i].second;
        s=vecOfNodes[s_i];
        L(s_i,s_i)-=1;
        reverseBool[c]=false;
        search = true;
        l= 0;
        while((l<Adj[s].size())&&(search)){
          if(Adj[s][l]==c){
            Adj[s].erase(Adj[s].begin()+l);
            search=false;
          }
          l++;
        }//end while
        for(count j=0;j<Adj[c].size();j++){
          w=Adj[c][j];
          w_i=reverse[w];
          if(reverseBool[w]&&(w!=s)){
            if(L(s_i,w_i)==0){
              L(s_i,w_i)=-1.;
              L(w_i,s_i)=-1.;
              L(s_i,s_i)+=1.;
              Adj[s].push_back (w);
              vecOfTransferredEdges.push_back(std::make_tuple(c,w,3));
              search = true;
              l= 0;
              while((l<Adj[w].size())&&(search)){
                if(Adj[w][l]==c){
                  Adj[w][l]=s;
                  search=false;
                }
                l++;
              }//end while
            }
            else{
              vecOfTransferredEdges.push_back(std::make_tuple(c,w,2));
              L(w_i,w_i)-=1.;
              search = true;
              l= 0;
              while((l<Adj[w].size())&&(search)){
                if(Adj[w][l]==c){
                  Adj[w].erase(Adj[w].begin()+l);
                  search=false;
                }
                l++;
              }//end while
            }//end else
          }//end if
        }//end for
        vecOfTransferredEdges.push_back(std::make_tuple(c,s,1));
        Adj[c].resize(0);
      }
      arma::uvec indices(vecOfNodes.size()-matchings.size());
      a=0;
      l=0;
      for(count i=0;i<vecOfNodes.size();i++){
        if(matchings[l].first==i){
          if(l<matchings.size()-1)
            l++;
        }
        else{
          indices(a)=i;
          a++;
        }
      }
      for(count i=0;i<matchings.size();i++){
        vecOfNodes.erase(vecOfNodes.begin()+matchings[i].first-i);
      }
      L=L.submat(indices, indices);

      ERDLevel Level(ID,vecOfNodes,vecOfTransferredEdges);
      LevelList.push_back(Level);
    }
} /* namespace NetworKit*/
