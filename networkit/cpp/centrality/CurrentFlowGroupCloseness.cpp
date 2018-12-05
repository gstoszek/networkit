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

   CurrentFlowGroupCloseness::CurrentFlowGroupCloseness(Graph G,const count k, const count CB,const double epsilon) : Centrality(G, true),k(k),CB(CB),epsilon(epsilon){
     if (G.isDirected()) throw std::runtime_error("Graph is directed graphs!");
     ConnectedComponents cc(G);
     cc.run();
     if (cc.getPartition().numberOfSubsets() > 1) throw std::runtime_error("Graph has more then one component!");
     if(k>=G.numberOfNodes()) throw std::runtime_error("Size of Group greater then number of nodes!");

     auto start = std::chrono::high_resolution_clock::now();
     auto end = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double> diff;

     numberOfCoarsedNodes=0;
     S.resize(k);
     n=G.numberOfNodes();

     end = std::chrono::high_resolution_clock::now();
     diff = end-start;
     std::cout << "Constructor finished in " << diff.count() << "(s)" << "\n\n";
           std::cout << "Number of Nodes=" << G.numberOfNodes() << "\n";
   }

   void CurrentFlowGroupCloseness::run() {
     bool coarse;
     count ID;
     count minDegree;
     auto start = std::chrono::high_resolution_clock::now();
     auto end = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double> diff;
     std::vector<std::tuple<count,count,count>> c_indices;
     arma::Mat<double> L;
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
         coarseGraph(c_indices,ID);
         ID++;
         minDegree=updateMinDegree();

         if(!(c_indices.size()>1)||!(G.numberOfNodes()>k)){
           coarse=false;
         }
        }
      }
      end = std::chrono::high_resolution_clock::now();
      diff = end-start;
      std::cout << "Coarsening " << diff.count() << "(s)" << "\n\n";
      std::cout << "Laplacian size afvoidter coarsening" << L.n_rows << "\n";
      start = std::chrono::high_resolution_clock::now();
      /*changed
      //ERD.computeFromLaplacian(vecOfNodes,L);
      */
      L=computePinvOfLaplacian();
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
      greedy(L);
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
    void CurrentFlowGroupCloseness::greedy(arma::Mat<double> L){
      count sampleSize;
      node s,v,w;
      double centrality,prevCFGCC,bestMarginalGain,distance;
      std::vector<bool> V;
      std::vector<node> vecOfNodes,vecOfSamples, vecOfPeriphs;
      std::vector<count> reverse;
      std::vector<double> mindst, dst, bst, zeroVec, marginalGain;

      vecOfSamples.resize(0);
      vecOfPeriphs.resize(0);
      vecOfNodes = G.nodes();
      reverse.resize(G.upperNodeIdBound());
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
      std::vector<node> vecOfNodes;
      node v;
      search=true;
      min=n;
      i=0;
      vecOfNodes=G.nodes();
      while((i<vecOfNodes.size())&&(search)){
        v=vecOfNodes[i];
        if((G.degree(v)<min) && (G.degree(v)>1)){
          min=G.degree(v);
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
        std::vector<node> vecOfNodes, vecOfCoarseNodes, vecOfNeighbours;
        std::vector<std::tuple<count,count,count>> indices;

        vecOfCoarseNodes.resize(0);
        vecOfFreeNodes.resize(G.upperNodeIdBound(),true);

        vecOfNodes=G.nodes();
        for(count i=0;i<vecOfNodes.size();i++){
          c=vecOfNodes[i];
          if(G.degree(c)==courseningDegree){
            vecOfCoarseNodes.push_back(i);
          }
        }
        if(Random){
          std::random_shuffle (vecOfCoarseNodes.begin(), vecOfCoarseNodes.end());
        }
        for(count i=0;i<vecOfCoarseNodes.size();i++){
          c=vecOfCoarseNodes[i];
          if(vecOfFreeNodes[c]){
            search=true;
            l=0;
            vecOfNeighbours=G.neighbors(c);
            while((search)&&(l<vecOfNeighbours.size())){
              s=vecOfNeighbours[l];
              if(vecOfFreeNodes[s]){
                l++;
                while((search)&&(l<vecOfNeighbours.size())){
                  w=vecOfNeighbours[l];
                  if(vecOfFreeNodes[w]){
                    vecOfFreeNodes[s]=false;
                    vecOfFreeNodes[c]=false;
                    vecOfFreeNodes[w]=false;
                    search=false;
                    indices.push_back(std::make_tuple(c,s,w));
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
    /**************************************************************************
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
        ERD.uncoarseTriangle(vecOfNodes,triangle);
        c=std::get<0>(triangle);
        s=std::get<1>(triangle);
        w=std::get<2>(triangle);
        edgeWeightcs=std::get<3>(triangle);
        edgeWeightcw=std::get<4>(triangle);
        edgeWeightsw=1./edgeWeightcs+1./edgeWeightsw;
        edgeWeightsw=1./edgeWeightsw;
        //vecOfNodes.push_back(c);
        //Adj[s].push_back (c);
        //Adj[c].push_back(s);
        if(std::get<4>(triangle)>0){
          //Adj[c].push_back (w);
          //Adj[w].push_back(c);
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
      count l;
      node c,s,w;
      double edgeWeightcs,edgeWeightss,edgeWeightsw;
      std::vector<node> vecOfNodes,vecOfSupernodes;
      std::vector<std::pair<node,node>> mapping;
      std::vector<std::tuple<node,node,node,double,double,double>> vecOfTriangles;

      vecOfNodes=G.nodes();
      vecOfSupernodes.resize(0);
      mapping.resize(0);
      vecOfTriangles.resize(0);
      for(count c_i=0;c_i<vecOfNodes.size();c_i++){
        c=vecOfNodes[c_i];
        if(G.degree(c)==1){
          s=G.randomNeighbor(c);
          if(std::find(vecOfSupernodes.begin(), vecOfSupernodes.end(), s) == vecOfSupernodes.end()){
            vecOfSupernodes.push_back(s);
            mapping.push_back (std::make_pair(s,c));
            vecOfPeripheralNodes[c]=true;
          }
          else{
            search = true;
            l= 0;
            while(search){
              if(vecOfSupernodes[l]==s){
                w=mapping[l].second;
                search=false;
              }
              l++;
            }//end while
            edgeWeightcs=G.weight(c,s);
            edgeWeightsw=G.weight(s,w);
            edgeWeightsw=1./edgeWeightcs+1./edgeWeightsw;
            edgeWeightsw=1./edgeWeightsw;
            G.setWeight(s,w,edgeWeightsw);
            G.removeNode(c);
            vecOfTriangles.push_back(std::make_tuple(c,s,w,edgeWeightcs,0.,edgeWeightsw));
          }//end else
        }
      }
      ERDLevel Level(1,vecOfTriangles);
      LevelList.push_back(Level);
    }
    /***************************************************************************/
    void CurrentFlowGroupCloseness::coarseGraph(std::vector<std::tuple<count,count,count>> matchings,count ID){
      bool search;
      count a,l,c_i,s_i,w_i;
      node c,s,w;
      double edgeWeightcs,edgeWeightcw,edgeWeightsw;
      std::vector<node> vecOfNodes;
      std::vector<std::tuple<node,node,node,double,double,double>> vecOfTriangles;

      vecOfNodes=G.nodes();
      vecOfTriangles.resize(0);

      for(count i=0;i<matchings.size();i++){
        search=false;
        c=std::get<0>(matchings[i]);
        s=std::get<1>(matchings[i]);
        w_i=std::get<2>(matchings[i]);
        w=vecOfNodes[w_i];
        edgeWeightcs=G.weight(c,s);
        edgeWeightcw=G.weight(c,w);
        edgeWeightsw=1./edgeWeightcs+1/edgeWeightcw;
        edgeWeightsw=1./edgeWeightsw;
        edgeWeightsw+=G.weight(s,w);
        G.setWeight(s,w,edgeWeightsw);
        vecOfTriangles.push_back(std::make_tuple(c,s,w,edgeWeightcs,edgeWeightcw,edgeWeightsw));
        G.removeNode(c);
        numberOfCoarsedNodes++;
      }
      /*
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
      */
      ERDLevel Level(ID,vecOfTriangles);
      LevelList.push_back(Level);
    }

    arma::Mat<double> CurrentFlowGroupCloseness::computePinvOfLaplacian(){
      node v;
      node w;
      count n;
      double factor;
      std::vector<node> vecOfNodes;
      vecOfNodes=G.nodes();
      n=vecOfNodes.size();
      arma::Mat<double> L(n,n);
      L.zeros();
      for(count i=0;i<n;i++){
        v=vecOfNodes[i];
        for(count j=i+1;j<n;j++){
          w=vecOfNodes[j];
          if(G.hasEdge(v,w)){
            L(i,j)=-G.weight(v,w);
            L(j,i)=-G.weight(v,w);
            L(i,i)+=G.weight(v,w);
            L(j,j)+=G.weight(v,w);
          }
        }
      }
      arma::Mat<double> J(n,n);
      factor=1./n;
      J.fill(factor);
      L= L+J;
      L=arma::inv_sympd(L);
      L= L-J;
      return L;
    }
} /* namespace NetworKit*/
