/*
 * CurrentFlowGroupCloseness.cpp
 *
 *      Author: gstoszek
 */
#define ARMA_DONT_PRINT_ERRORS
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

   CurrentFlowGroupCloseness::CurrentFlowGroupCloseness(Graph& G,const count k, const count CB,const double epsilon, const bool doInvert)
    : G(G),k(k),CB(CB),epsilon(epsilon),doInvert(doInvert){
     if (G.isDirected()) throw std::runtime_error("Graph is directed!");
     ConnectedComponents cc(G);
     cc.run();
     if (cc.getPartition().numberOfSubsets() > 1) throw std::runtime_error("Graph has more then one component!");
     if(k>=G.numberOfNodes()) throw std::runtime_error("Size of Group greater then number of nodes!");
     auto start = std::chrono::high_resolution_clock::now();
     auto end = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double> diff;

     S.resize(k);
     n=G.numberOfNodes();

     end = std::chrono::high_resolution_clock::now();
     diff = end-start;
     std::cout << "Constructor finished in " << diff.count() << "(s)" << "\n";
     std::cout << "Number of Nodes=" << G.numberOfNodes() << "\n";
   }

   void CurrentFlowGroupCloseness::run() {
     bool coarse;
     count ID, minDegree;
     auto start = std::chrono::high_resolution_clock::now();
     auto end = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double> diff;
     std::vector<node> vecOfChosenNodes;

     ID=0;
     if(CB>0){
       ID++;
       minDegree=1;
       mergePeripheralNodes();
       ID++;
       minDegree=updateMinDegree();
       coarse=true;
       while((minDegree<CB)&&(coarse)){
         vecOfChosenNodes=coarsingIndices(minDegree, false);
         coarseGraph(vecOfChosenNodes,ID);
         ID++;
         minDegree=updateMinDegree();
         if(!(vecOfChosenNodes.size()>1)||!(G.numberOfNodes()>k)){
           coarse=false;
         }
        }
      }
      end = std::chrono::high_resolution_clock::now();
      diff = end-start;
      std::cout << "Coarsening of G finished in: " << diff.count() << "(s)" << "\n";
      std::cout << "Number of nodes after coarsening: " << G.numberOfNodes() << "\n";

      /*
      start = std::chrono::high_resolution_clock::now();
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
      end = std::chrono::high_resolution_clock::now();
      diff = end-start;
      std::cout << "Uncoarsening in: " << diff.count() << "(s)" << "\n";
      */
      std::cout << "Starting Greedy-Algorithm\n";
      start = std::chrono::high_resolution_clock::now();
      if(doInvert)
        greedy();
      else{
        greedyLAMG();
      }
      end = std::chrono::high_resolution_clock::now();
      diff = end-start;
      std::cout << "Greedy in: " << diff.count() << "(s)" << "\n\n";
    }
    std::vector<node> CurrentFlowGroupCloseness::getNodesofGroup(){
      return S;
    }
    double CurrentFlowGroupCloseness::getCFGCC() {
      return CFGCC;
    }
    void CurrentFlowGroupCloseness::greedy(){
      count v_i,w_i;
      node s,v,w;
      double centrality,prevCFGCC,bestMarginalGain,distance;
      std::vector<bool> V,W;
      std::vector<node> vecOfNodes,vecOfSamples, vecOfPeriphs,cutNodes;
      std::vector<count> reverse;
      std::vector<double> mindst, dst, bst, minApprox, bstApprox, marginalGain;

      arma::Mat<double> Pinv;
      Pinv=computePinvOfLaplacian();

      CFGCC=n*n*n;
      prevCFGCC=CFGCC;
      V.resize(G.numberOfNodes(),true);
      W.resize(G.numberOfNodes(),true);
      mindst.resize(G.numberOfNodes(),n*n);
      marginalGain.resize(G.numberOfNodes(),CFGCC);

      vecOfSamples.resize(0);
      vecOfPeriphs.resize(0);
      vecOfNodes = G.nodes();
      reverse.resize(G.upperNodeIdBound());
      for(count i=0;i<vecOfNodes.size();i++){
        reverse[vecOfNodes[i]]=i;
      }
      for(count i=0;i<vecOfNodes.size();i++){
        v=vecOfNodes[i];
        if(std::find(vecOfPeripheralNodes.begin(), vecOfPeripheralNodes.end(), v) != vecOfPeripheralNodes.end())
          vecOfSamples.push_back(v);
        else{
          vecOfPeriphs.push_back(v);
        }
      }
      for(count i=0;i<k;i++){
        std::random_shuffle (vecOfSamples.begin(), vecOfSamples.end());
        bestMarginalGain=0.;
        for (count j=0; j<vecOfSamples.size();j++) {
          v=vecOfSamples[j];
          v_i=reverse[v];
          if(V[v_i] && (bestMarginalGain<marginalGain[v_i])){
            dst=mindst;
            centrality = 0.;
            cutNodes.resize(0);
            for (count l = 0; l < vecOfSamples.size(); l++) {
              w=vecOfSamples[l];
              w_i=reverse[w];
              if(W[w_i]){
                distance=Pinv(v_i,v_i)+Pinv(w_i,w_i)-2*Pinv(v_i,w_i);
                if (distance< mindst[w_i]){
                  dst[w_i]=distance;
                  if(distance<bst[w_i])
                    cutNodes.push_back(w_i);
                }
              }
              centrality +=dst[w_i];
            }
            for (count l = 0 ; l < vecOfPeriphs.size(); l++) {
              w=vecOfPeriphs[l];
              w_i=reverse[w];
              if(W[w_i]){
                distance=Pinv(v_i,v_i)+Pinv(w_i,w_i)-2*Pinv(v_i,w_i);
                if (distance< mindst[w_i]){
                  dst[w_i]=distance;
                  if(distance<bst[w_i])
                    cutNodes.push_back(w_i);
                }
                centrality += dst[w_i];
              }
            }
            marginalGain[v_i]=prevCFGCC-centrality;
            if (centrality < CFGCC) {
              CFGCC = centrality;
              bst=dst;
              s = v;
              bestMarginalGain=marginalGain[v_i];
            }
          }
        }
        S[i]=s;
        V[s]=false;
        mindst=bst;
        prevCFGCC=CFGCC;
      }
     CFGCC = (double)(n)/CFGCC;
    }
    void CurrentFlowGroupCloseness::greedyLAMG(){
      count v_i,w_i;
      node s,v,w;
      double centrality,prevCFGCC,bestMarginalGain,distance;
      std::vector<bool> V;
      std::vector<node> vecOfNodes,vecOfSamples, vecOfPeriphs;
      std::vector<count> reverse;
      std::vector<double> mindst, dst, bst, minApprox, bstApprox, marginalGain;
      std::cout<<"***************Setup***************\n";
      Lamg<CSRMatrix> lamg;
      CSRMatrix matrix = CSRMatrix::laplacianMatrix(G);
      lamg.setupConnected(matrix);
      std::cout<<"**************Setup***************\n";
      CFGCC=n*n*n;
      prevCFGCC=CFGCC;
      mindst.resize(G.numberOfNodes(),n*n);
      marginalGain.resize(G.numberOfNodes(),CFGCC);
      V.resize(G.numberOfNodes(),true);
      vecOfSamples.resize(0);
      vecOfPeriphs.resize(0);
      vecOfNodes = G.nodes();
      reverse.resize(G.upperNodeIdBound());
      for(count i=0;i<vecOfNodes.size();i++){
        reverse[vecOfNodes[i]]=i;
      }
      for(count i=0;i<vecOfNodes.size();i++){
        v=vecOfNodes[i];
        if(std::find(vecOfPeripheralNodes.begin(), vecOfPeripheralNodes.end(), v) != vecOfPeripheralNodes.end())
          vecOfSamples.push_back(v);
        else{
          vecOfPeriphs.push_back(v);
        }
      }
      Vector result(vecOfNodes.size());
      Vector rhs(vecOfNodes.size(), 0.);
      Vector zeroVector(vecOfNodes.size(), 0.);
      std::cout<<"********************1***************\n";
      for(count i=0;i<k;i++){
        std::random_shuffle (vecOfSamples.begin(), vecOfSamples.end());
        bestMarginalGain=0.;
        for (count j=0; j<vecOfSamples.size();j++) {
          v=vecOfSamples[j];
          v_i=reverse[v];
          if(V[v_i] && (bestMarginalGain<marginalGain[v_i])){
            rhs[v_i]=1.;
            dst=mindst;
            centrality = 0.;
            for (count l = 0; l < vecOfSamples.size(); l++) {
              if(v!=w){
                w=vecOfSamples[l];
                w_i=reverse[w];
                rhs[w_i]=-1.;
                result=zeroVector;
                lamg.solve(rhs, result);
                distance=fabs(result[v_i]-result[w_i]);
                if (distance< mindst[w_i])
                  dst[w_i]=distance;
                centrality +=dst[w_i];
                rhs[w_i]=0.;
              }
            }
            for (count l = 0 ; l < vecOfPeriphs.size(); l++) {
              w=vecOfPeriphs[l];
              w_i=reverse[w];
              rhs[w_i]=-1.;
              result=zeroVector;
              lamg.solve(rhs, result);
              distance=fabs(result[v_i]-result[w_i]);
              if (distance< mindst[w_i])
                dst[w_i]=distance;
              centrality += dst[w_i];
              rhs[w_i]=0.;
            }
            marginalGain[v_i]=prevCFGCC-centrality;
            if (centrality < CFGCC) {
              CFGCC = centrality;
              bst=dst;
              s = v;
              bestMarginalGain=marginalGain[v_i];
            }
            rhs[v_i]=0.;
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

    std::vector<node> CurrentFlowGroupCloseness::coarsingIndices(count courseningDegree, bool Random){
        bool search;
        count l;
        node c,w;
        std::vector<node> vecOfNodes, vecOfChosenNodes, vecOfNeighbors,vecOfOccupiedNodes,vecOfCandidates;

        vecOfOccupiedNodes.resize(0);
        vecOfCandidates.resize(0);
        vecOfNodes=G.nodes();
        for(count i=0;i<vecOfNodes.size();i++){
          c=vecOfNodes[i];
          if(G.degree(c)==courseningDegree){
            vecOfCandidates.push_back(c);
          }
        }
        if(Random){
          std::random_shuffle (vecOfCandidates.begin(), vecOfCandidates.end());
        }
        for(count i=0;i<vecOfCandidates.size();i++){
          c=vecOfCandidates[i];
          if(std::find(vecOfOccupiedNodes.begin(), vecOfOccupiedNodes.end(), c) == vecOfOccupiedNodes.end()){
            vecOfNeighbors=G.neighbors(c);
            search=true;
            l=0;
            if((search)&&(l<vecOfNeighbors.size())){
              w=vecOfNeighbors[l];
              if(std::find(vecOfOccupiedNodes.begin(), vecOfOccupiedNodes.end(), w) != vecOfOccupiedNodes.end()){
                search=false;
              }
              else{
                l++;
              }
            }
            if(search){
              for(count j=0;j<vecOfNeighbors.size();j++){
                w=vecOfNeighbors[j];
                vecOfOccupiedNodes.push_back(w);
              }
              vecOfChosenNodes.push_back(c);
            }
          }
        }
        return vecOfChosenNodes;
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
      double edgeWeightcs,edgeWeightsw;
      std::vector<node> vecOfNodes,vecOfSupernodes;
      std::vector<double> vecOfWeights;
      std::vector<std::pair<node,node>> mapping;

      vecOfSupernodes.resize(0);
      mapping.resize(0);
      vecOfNodes=G.nodes();
      vecOfWeights.resize(1);
      for(count i=0;i<vecOfNodes.size();i++){
        c=vecOfNodes[i];
        if(G.degree(c)==1){
          s=G.randomNeighbor(c);
          if(std::find(vecOfSupernodes.begin(), vecOfSupernodes.end(), s) == vecOfSupernodes.end()){
            vecOfSupernodes.push_back(s);
            mapping.push_back (std::make_pair(s,c));
            vecOfPeripheralNodes.push_back(c);
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
            edgeWeightsw=edgeWeightcs+edgeWeightsw;
            G.setWeight(s,w,edgeWeightsw);
            vecOfWeights[0]=edgeWeightcs;
            ERDLevel coarsed(c,G.neighbors(c),vecOfWeights);
            coarsedNodes.push_back(coarsed);
            G.removeNode(c);
          }//end else
        }
      }
    }
    /***************************************************************************/
    void CurrentFlowGroupCloseness::coarseGraph(std::vector<node> vecOfChosenNodes,count degree){
      node c,s,w;
      double edgeWeightcs,edgeWeightcw,edgeWeightsw;
      std::vector<node> vecOfNeighbors;
      std::vector<double> vecOfWeights;

      vecOfWeights.resize(degree);
      /*Path Coarsenening*/
      if(degree==2){
        for(count i=0;i<vecOfChosenNodes.size();i++){
          c=vecOfChosenNodes[i];
          vecOfNeighbors=G.neighbors(c);
          s=vecOfNeighbors[0];
          w=vecOfNeighbors[1];
          edgeWeightcs=G.weight(c,s);
          edgeWeightcw=G.weight(c,w);
          edgeWeightsw=edgeWeightcs+edgeWeightcw;
          /*Path Merge*/
          if(G.weight(s,w)!=0){
            edgeWeightsw=1./edgeWeightsw;
            edgeWeightsw+=1./G.weight(s,w);
            edgeWeightsw=1./edgeWeightsw;
          }
          G.setWeight(s,w,edgeWeightsw);
          for(count j=0;j<vecOfNeighbors.size();j++){
            vecOfWeights[j]=G.weight(c,vecOfNeighbors[j]);
          }
          ERDLevel coarsed(c,vecOfNeighbors,vecOfWeights);
          G.removeNode(c);
        }
      }
      /*Star*/
      else{
        for(count i=0;i<vecOfChosenNodes.size();i++){
          c=vecOfChosenNodes[i];
          computeStarCliqueWeights(c);
          vecOfNeighbors=G.neighbors(c);
          for(count j=0;j<vecOfNeighbors.size();j++){
            vecOfWeights[j]=G.weight(c,vecOfNeighbors[j]);
          }
          ERDLevel coarsed(c,vecOfNeighbors,vecOfWeights);
          G.removeNode(c);
        }
      }
    }

    void CurrentFlowGroupCloseness::computeStarCliqueWeights(node c){
      node v,w,x,y;
      double S,S2,weight;
      std::vector<node> vecOfNeighbors;

      S=0.;
      vecOfNeighbors=G.neighbors(c);
      for(count i=0;i<vecOfNeighbors.size();i++){
        v=vecOfNeighbors[i];
        S2=0.;
        for(count j=0;j<vecOfNeighbors.size();j++){
          w=vecOfNeighbors[j];
          if(w!=v){
            S2*=G.weight(v,w);
          }
        }
        S+=S2;
      }
      for(count i=0;i<vecOfNeighbors.size();i++){
        v=vecOfNeighbors[i];
        for(count j=0;j<vecOfNeighbors.size();j++){
          w=vecOfNeighbors[j];
          if(v!=w){
            weight=0.;
            for(count k=0;k<vecOfNeighbors.size();k++){
              x=vecOfNeighbors[k];
              if((x!=v)&&(x!=w)){
                for(count l=0;l<vecOfNeighbors.size();l++){
                  y=vecOfNeighbors[l];
                  if((y!=v)&&(y!=w)&&(y!=x)){
                    weight*=G.weight(x,y);
                  }//if y
                }//for y
              }//if x
            }//for x
            weight=S/weight;
            if(G.weight(v,w)!=0){
              weight=1./weight;
              weight+=1./G.weight(v,w);
              weight=1./weight;
            }
            G.setWeight(x,y,weight);
          }// if w
        }// for w
      }// for v
    }
    arma::Mat<double> CurrentFlowGroupCloseness::computePinvOfLaplacian(){
      node v,w;
      count w_i;
      double factor, weight;
      std::vector<node> vecOfNodes, vecOfNeighbours;
      std::vector<count> reverse;
      vecOfNodes=G.nodes();
      arma::Mat<double> L(G.numberOfNodes(),G.numberOfNodes());
      L.zeros();
      reverse.resize(G.upperNodeIdBound());
      for(count i=0;i<vecOfNodes.size();i++){
        reverse[vecOfNodes[i]]=i;
      }
      for(count i=0;i<G.numberOfNodes();i++){
        v=vecOfNodes[i];
        vecOfNeighbours=G.neighbors(v);
        for(count j=0;j<vecOfNeighbours.size();j++){
          w=vecOfNeighbours[j];
          w_i=reverse[w];
          weight=G.weight(v,w);
          L(i,w_i)=-weight;
          L(i,i)+=weight;
        }
      }
      arma::Mat<double> J(L.n_rows,L.n_rows);
      factor=1./(double)(L.n_rows);
      J.fill(factor);
      L= L-J;
      L=arma::inv(L);
      L= L+J;
      return L;
    }
} /* namespace NetworKit*/
