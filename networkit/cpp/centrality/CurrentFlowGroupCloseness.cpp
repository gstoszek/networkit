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

     auto start = std::chrono::high_resolution_clock::now();
     auto end = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double> diff;

     node v;
     node w;
     double edgeWeight;
     numberOfCoarsedNodes=0;
     S.resize(k);
     n=G.numberOfNodes();
     Adj.resize(n);
     L.set_size(n,n);
     L.zeros();
     std::vector<node> vecOfNodesMapping;
     vecOfPeripheralNodes.resize(n,false);
     vecOfNodes.resize(n);
     vecOfNodesMapping=G.nodes();
     //ERD.M=L;
     for(count i=0;i<L.n_rows;i++){
       vecOfNodes[i]=i;
       v=vecOfNodesMapping[i];
       for(count j=i+1;j<L.n_rows;j++){
         w=vecOfNodesMapping[j];
         if(G.hasEdge(v,w)){
           edgeWeight=1./G.weight(v,w);
           L(i,j)=-1*edgeWeight;
           L(j,i)=L(i,j);
           L(i,i)+=edgeWeight;
           L(j,j)+=edgeWeight;
           Adj[i].push_back(j);
           Adj[j].push_back(i);
         }
       }
     }
     end = std::chrono::high_resolution_clock::now();
     diff = end-start;
     std::cout << "Constructor finished in " << diff.count() << "(s)" << "\n\n";
           std::cout << "Laplacian" << L.n_rows << "\n";
   }

   void CurrentFlowGroupCloseness::run() {
     bool coarse;
     count ID;
     count minDegree;
     auto start = std::chrono::high_resolution_clock::now();
     auto end = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double> diff;
     std::vector<std::tuple<count,count,count>> c_indices;;
     ID=0;
     if(CB>0){
       ID++;
       minDegree=1;
       mergePeripheralNodes();
       ID++;
       minDegree=updateMinDegree();
       coarse=true;
       while((minDegree<CB)&&(coarse)){
         c_indices=coarsingIndices(minDegree, false);
         coarseLaplacian(c_indices,ID);
         ID++;
         minDegree=updateMinDegree();
         if(!(c_indices.size()>1)||!(vecOfNodes.size()>k)){
           coarse=false;
         }
        }
      }
      end = std::chrono::high_resolution_clock::now();
      diff = end-start;
      std::cout << "Coarsening " << diff.count() << "(s)" << "\n\n";
      std::cout << "Laplacian size after coarsening" << L.n_rows << "\n";
      start = std::chrono::high_resolution_clock::now();
      /*changed
      //ERD.computeFromLaplacian(vecOfNodes,L);
      */
      computePinvOfLaplacian();
      end = std::chrono::high_resolution_clock::now();
      diff = end-start;
      std::cout << "Initial EffectiveResistanceDistanceMatrix finished in " << diff.count() << "(s)" << "\n\n";
      //greedy();
      //std::cout<<"Level:" <<ID<<" with the value: " << CFGCC << "\n";
      start = std::chrono::high_resolution_clock::now();
      /*
      if(CB>1){
        ID=LevelList[LevelList.size()-1].getID();
        std::cout<<"ID="<<ID<<"\n";
        while(ID>1){
          uncoarseEfffectiveResistanceDistanceMatrix(ID);
          //greedy();
          //std::cout<<"Level:" <<ID<<" with the value: " << CFGCC << "\n";
          ID--;
        }
      }
      */
      end = std::chrono::high_resolution_clock::now();
      diff = end-start;
      std::cout << "Uncoarsening in " << diff.count() << "(s)" << "\n\n";
      start = std::chrono::high_resolution_clock::now();
      greedy();
      end = std::chrono::high_resolution_clock::now();
      diff = end-start;
      std::cout << "Greedy in " << diff.count() << "(s)" << "\n\n";
    }
    std::vector<node> CurrentFlowGroupCloseness::getNodesofGroup(){
      return S;
    }
    double CurrentFlowGroupCloseness::getCFGCC() {
      return CFGCC;
    }
    void CurrentFlowGroupCloseness::greedy(){
      count sampleSize;
      node s,v,w;
      double centrality,prevCFGCC,bestMarginalGain,distance;
      std::vector<bool> V;
      std::vector<node> vecOfSamples, vecOfPeriphs;
      std::vector<count> reverse;
      std::vector<double> mindst, dst, bst, zeroVec, marginalGain;

      vecOfSamples.resize(0);
      vecOfPeriphs.resize(0);

      reverse.resize(G.numberOfNodes());
      for(count i=0;i<vecOfNodes.size();i++){
        reverse[vecOfNodes[i]]=i;
      }

      for(count i=0;i<vecOfNodes.size();i++){
        v=vecOfNodes[i];
        if(vecOfPeripheralNodes[v])
          vecOfPeriphs.push_back(v);
        else{
          vecOfSamples.push_back(v);
        }
      }
      sampleSize=(count)(log(vecOfSamples.size())/(2*epsilon*epsilon));
      if(sampleSize>vecOfSamples.size()){
        sampleSize=vecOfSamples.size();
      }
      CFGCC=n*n*n;
      prevCFGCC=CFGCC;
      V.resize(G.numberOfNodes(),true);
      mindst.resize(G.numberOfNodes(),n*n);
      zeroVec.resize(G.numberOfNodes(),0.);
      marginalGain.resize(G.numberOfNodes(),CFGCC);
      for(count i=0;i<k;i++){
        std::random_shuffle (vecOfSamples.begin(), vecOfSamples.end());
        bestMarginalGain=0.;
        for (count j=0; j<vecOfSamples.size();j++) {
          v=vecOfSamples[j];
          if(V[v] && (bestMarginalGain<marginalGain[v])){
            dst=mindst;
            centrality = 0.;
            for (count l = 0; l < sampleSize; l++) {
              w=vecOfSamples[l];
              distance=L(reverse[v],reverse[v])+L(reverse[w],reverse[w])-2*L(reverse[v],reverse[w]);
              if (distance< mindst[w]){
                dst[w]=distance;
              }
              centrality +=dst[w];
            }
            centrality*=((double)(vecOfSamples.size()+numberOfCoarsedNodes)/(double)(sampleSize));
            for (count l = 0 ; l < vecOfPeriphs.size(); l++) {
              w=vecOfPeriphs[l];
              distance=L(reverse[v],reverse[v])+L(reverse[w],reverse[w])-2*L(reverse[v],reverse[w]);
              if (distance< mindst[w]){
                dst[w]=distance;
              }
              centrality += dst[w];
            }
            marginalGain[v]=prevCFGCC-centrality;
            if (centrality < CFGCC) {
              CFGCC = centrality;
              bst=dst;
              s = v;
              bestMarginalGain=marginalGain[v];
            }
            //std::cout<<"c("<<v<<")"<<centrality<<"\n";
          }
        }
        S[i]=s;
        V[s]=false;
        mindst=bst;
        prevCFGCC=CFGCC;
      }
     CFGCC = (double)(n)/CFGCC;
    }
    count CurrentFlowGroupCloseness::updateMinDegree(){
      bool search;
      count min;
      count i;
      node v;
      search=true;
      min=n;
      i=0;

      while((i<L.n_rows)&&(search)){
        v=vecOfNodes[i];
        if((Adj[v].size()<min) && (Adj[v].size()>1)){
          min=Adj[v].size();
          if(min==2){
            search=false;
          }
        }
        i++;
      }
      return min;
    }
    std::vector<std::tuple<count,count,count>> CurrentFlowGroupCloseness::coarsingIndices(count courseningDegree, bool Random){
        bool search;
        count c_i,s_i,w_i,l;
        node c,s,w;
        std::vector<bool> vecOfFreeNodes;
        std::vector<count> vecOfCoarseNodes;
        std::vector<std::tuple<count,count,count>> indices;
        std::vector<count> reverse;

        vecOfCoarseNodes.resize(0);
        vecOfFreeNodes.resize(G.numberOfNodes(),true);
        reverse.resize(G.numberOfNodes());
        for(count i=0;i<vecOfNodes.size();i++){
          reverse[vecOfNodes[i]]=i;
        }
        for(count i=0;i<L.n_rows;i++){
          c=vecOfNodes[i];
          if(Adj[c].size()==courseningDegree){
            vecOfCoarseNodes.push_back(i);
          }
        }
        if(Random){
          std::random_shuffle (vecOfCoarseNodes.begin(), vecOfCoarseNodes.end());
        }
        for(count i=0;i<vecOfCoarseNodes.size();i++){
          c_i=vecOfCoarseNodes[i];
          c=vecOfNodes[c_i];
          if(vecOfFreeNodes[c]){
            search=true;
            l=0;
            while((search)&&(l<Adj[c].size())){
              s=Adj[c][l];
              if(vecOfFreeNodes[s]){
                l++;
                while((search)&&(l<Adj[c].size())){
                  w=Adj[c][l];
                  if(vecOfFreeNodes[w]){
                    vecOfFreeNodes[s]=false;
                    vecOfFreeNodes[c]=false;
                    vecOfFreeNodes[w]=false;
                    s_i=reverse[s];
                    w_i=reverse[w];
                    search=false;
                    indices.push_back(std::make_tuple(c_i,s_i,w_i));
                  }
                  else{
                    l++;
                  }
                }
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
      count l,j;
      node c,s,w;
      double edgeWeightcs,edgeWeightcw,edgeWeightsw;
      std::vector<std::tuple<node,node,count,double,double,double>> vecOfTriangles;
      std::tuple<node,node,count,double,double,double> triangle;

      vecOfTriangles=LevelList[ID-1].getVecOfTriangles();
      for(count i=0;i<vecOfTriangles.size();i++){
        j=vecOfTriangles.size()-1-i;
        triangle=vecOfTriangles[j];
        /**
        ERD.uncoarseTriangle(vecOfNodes,triangle);
        **/
        c=std::get<0>(triangle);
        s=std::get<1>(triangle);
        w=std::get<2>(triangle);
        edgeWeightcs=std::get<3>(triangle);
        edgeWeightcw=std::get<4>(triangle);
        edgeWeightsw=1./edgeWeightcs+1./edgeWeightsw;
        edgeWeightsw=1./edgeWeightsw;
        vecOfNodes.push_back(c);
        Adj[s].push_back (c);
        Adj[c].push_back(s);
        if(std::get<4>(triangle)>0){
          Adj[c].push_back (w);
          Adj[w].push_back(c);
        }
        if(!(std::get<5>(triangle)==edgeWeightsw)){
          search = true;
          l= 0;
          while((l<Adj[s].size())&&(search)){
            if(Adj[s][l]==w){
              Adj[s].erase(Adj[s].begin()+l);
              search=false;
            }
            else{
              l++;
            }
          }//end while
          search = true;
          l= 0;
          while((l<Adj[w].size())&&(search)){
            if(Adj[w][l]==s){
                Adj[w].erase(Adj[w].begin()+l);
                search=false;
              }
          else{
            l++;
          }
          }//end while
        }//end if
        numberOfCoarsedNodes--;
      }//end for
    }
    /***************************************************************************/
    void CurrentFlowGroupCloseness::mergePeripheralNodes(){
      bool search;
      count a,l,c_i,s_i,w_i;
      node c,s,w;
      double edgeWeightcs,edgeWeightsw;
      std::vector<bool> vecOfFreeNeighbours;
      std::vector<node> mapping;
      std::vector<node> reverse;
      std::vector<std::tuple<node,node,node,double,double,double>> vecOfTriangles;

      vecOfFreeNeighbours.resize(G.numberOfNodes(),true);
      mapping.resize(G.numberOfNodes());
      reverse.resize(G.numberOfNodes());
      for(count i=0;i<vecOfNodes.size();i++){
        reverse[vecOfNodes[i]]=i;
      }
      vecOfTriangles.resize(0);
      for(count c_i=0;c_i<L.n_rows;c_i++){
        c=vecOfNodes[c_i];
        if(Adj[c].size()==1){
          s=Adj[c][0];
          if(vecOfFreeNeighbours[s]){
            vecOfFreeNeighbours[s]=false;
            mapping[s]=c;
            vecOfPeripheralNodes[c]=true;
          }
          else{
            w=mapping[s];
            w_i=reverse[w];
            s_i=reverse[s];
            edgeWeightcs=-L(c_i,s_i);
            edgeWeightsw=-L(s_i,w_i);
            L(s_i,s_i)-=(edgeWeightsw+edgeWeightcs);
            edgeWeightsw=1./edgeWeightcs+1./edgeWeightsw;
            edgeWeightsw=1./edgeWeightsw;
            L(s_i,w_i)=-edgeWeightsw;
            L(w_i,s_i)=L(s_i,w_i);
            L(s_i,s_i)+=edgeWeightsw;
            L(w_i,w_i)=edgeWeightsw;
            search = true;
            l= 0;
            while((l<Adj[s].size())&&(search)){
              if(Adj[s][l]==c){
                Adj[s].erase(Adj[s].begin()+l);
                search=false;
              }
              l++;
            }//end while
            Adj[c].resize(0);
            vecOfTriangles.push_back(std::make_tuple(c,s,w,edgeWeightcs,0.,edgeWeightsw));
          }//end else
        }
      }
      if(vecOfTriangles.size()>0){
        arma::uvec indices(vecOfNodes.size()-vecOfTriangles.size());
        a=0;
        l=0;
        for(count i=0;i<vecOfNodes.size();i++){
          c=std::get<0>(vecOfTriangles[l]);
          if(reverse[c]==i){
            if(l<vecOfTriangles.size()-1)
            l++;
          }
          else{
            indices(a)=i;
            a++;
          }
        }
        for(count i=0;i<vecOfTriangles.size();i++){
          c=std::get<0>(vecOfTriangles[i]);
          c_i=reverse[c];
          vecOfNodes.erase(vecOfNodes.begin()+c_i-i);
        }
        L=L.submat(indices, indices);
      }
      ERDLevel Level(1,vecOfTriangles);
      LevelList.push_back(Level);
    }
    /***************************************************************************/
    void CurrentFlowGroupCloseness::coarseLaplacian(std::vector<std::tuple<count,count,count>> matchings,count ID){
      bool search;
      count a,l,c_i,s_i,w_i;
      node c,s,w;
      double edgeWeightcs,edgeWeightcw,edgeWeightsw;
      std::vector<count> reverse;
      std::vector<std::tuple<node,node,node,double,double,double>> vecOfTriangles;

      reverse.resize(G.numberOfNodes());
      vecOfTriangles.resize(0);
      for(count i=0;i<vecOfNodes.size();i++){
        reverse[vecOfNodes[i]]=i;
      }
      for(count i=0;i<matchings.size();i++){
        search=false;
        c_i=std::get<0>(matchings[i]);
        c=vecOfNodes[c_i];
        s_i=std::get<1>(matchings[i]);
        s=vecOfNodes[s_i];
        w_i=std::get<2>(matchings[i]);
        w=vecOfNodes[w_i];
        edgeWeightcs=-L(c_i,s_i);
        edgeWeightcw=-L(c_i,w_i);
        edgeWeightsw=-L(s_i,w_i);
        if(edgeWeightsw==0){
          search=true;
        }
        L(s_i,s_i)-=(edgeWeightcs+edgeWeightsw);
        L(w_i,w_i)-=(edgeWeightcw+edgeWeightsw);
        edgeWeightsw=1./edgeWeightcs+1/edgeWeightcw;
        edgeWeightsw=1./edgeWeightsw;
        L(s_i,w_i)-=edgeWeightsw;
        L(w_i,s_i)=L(s_i,w_i);
        edgeWeightsw=-L(s_i,w_i);
        L(s_i,s_i)+=edgeWeightsw;
        L(w_i,w_i)+=edgeWeightsw;
        vecOfTriangles.push_back(std::make_tuple(c,s,w,edgeWeightcs,edgeWeightcw,edgeWeightsw));
        if(search){
          Adj[s].push_back(w);
          Adj[w].push_back(s);
        }
        search=true;
        l= 0;
        while((l<Adj[s].size())&&(search)){
          if(Adj[s][l]==c){
            Adj[s].erase(Adj[s].begin()+l);
            search=false;
          }
          l++;
        }//end while
        search = true;
        l= 0;
        while((l<Adj[w].size())&&(search)){
          if(Adj[w][l]==c){
            Adj[w].erase(Adj[w].begin()+l);
            search=false;
          }
          l++;
        }//end while
        Adj[c].resize(0);
        numberOfCoarsedNodes++;
      }
      arma::uvec indices(vecOfNodes.size()-matchings.size());
      a=0;
      l=0;
      for(count i=0;i<vecOfNodes.size();i++){
        if(std::get<0>(matchings[l])==i){
          if(l<matchings.size()-1)
            l++;
        }
        else{
          indices(a)=i;
          a++;
        }
      }
      for(count i=0;i<matchings.size();i++){
        vecOfNodes.erase(vecOfNodes.begin()+std::get<0>(matchings[i])-i);
      }
      L=L.submat(indices, indices);
      ERDLevel Level(ID,vecOfTriangles);
      LevelList.push_back(Level);
    }

    void CurrentFlowGroupCloseness::computePinvOfLaplacian(){
      node v;
      node w;
      count n;
      double factor;
      n=L.n_rows;
      arma::Mat<double> J(n,n);
      factor=1./n;
      J.fill(factor);
      L= L+J;
      L=arma::inv_sympd(L);
      L= L-J;
    }
} /* namespace NetworKit*/
