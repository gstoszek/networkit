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
     limit=G.upperNodeIdBound();
     coarsedNodes.resize(0);
     coarsedNeighbors.resize(0);
     coarsedWeights.resize(0);

     end = std::chrono::high_resolution_clock::now();
     diff = end-start;
     std::cout << "Constructor finished in " << diff.count() << "(s)" << "\n";
     std::cout << "Number of Nodes=" << G.numberOfNodes() << "\n";
   }

   void CurrentFlowGroupCloseness::run() {
     bool coarse;
     count ID, degree;
     auto start = std::chrono::high_resolution_clock::now();
     auto end = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double> diff;
     std::vector<node> vecOfChosenNodes;

     ID=0;
     if(CB>0){
       ID++;
       degree=1;
       mergePeripheralNodes();
       ID++;
       degree=updateMinDegree();
       coarse=true;
       while((degree<CB)&&(coarse)){
         vecOfChosenNodes=coarsingIndices(degree, false);
         coarseGraph(vecOfChosenNodes,degree);
         ID++;
         degree=updateMinDegree();
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
      double d,dApprox,prevD,bestMarginalGain,distance;
      std::vector<bool> V,E;
      std::vector<node> vecOfNodes;
      std::vector<count> reverse;
      std::vector<double> Sdst, dst, tmpdst, SdstApprox, dstApprox, tmpApprox, marginalGain;

      arma::Mat<double> Pinv;
      Pinv=computePinvOfLaplacian();
      CFGCC=n*n*n;
      prevD=CFGCC;
      V.resize(limit,false);
      E=V;
      Sdst.resize(G.numberOfNodes(),n*n);
      SdstApprox.resize(coarsedNodes.size(),n*n);
      marginalGain.resize(G.numberOfNodes(),CFGCC);

      vecOfNodes = G.nodes();
      reverse.resize(limit);

      for(count i=0;i<coarsedNodes.size();i++){
        v=coarsedNodes[i];
        reverse[v]=i;
      }
      for(count i=0;i<vecOfNodes.size();i++){
        v=vecOfNodes[i];
        reverse[v]=i;
        E[v]=true;
        V[v]=true;
      }
      for(count i=0;i<vecOfPeripheralNodes.size();i++){
        V[vecOfPeripheralNodes[i]]=false;
      }
      for(count i=0;i<k;i++){
        bestMarginalGain=0.;
        for (count j=0; j<vecOfNodes.size();j++) {
          v=vecOfNodes[j];
          v_i=reverse[v];
          if(V[v] && (bestMarginalGain<marginalGain[v_i])){
            dst=Sdst;
            d = 0.;
            for (count l = 0; l < vecOfNodes.size(); l++) {
              w=vecOfNodes[l];
              w_i=reverse[w];
              distance=Pinv(v_i,v_i)+Pinv(w_i,w_i)-2*Pinv(v_i,w_i);
              dst[w_i]=std::min(distance,dst[w_i]);
              d +=dst[w_i];
            }
            dstApprox=SdstApprox;
            dApprox=computeApproxDistances(E, reverse, dst, &dstApprox);
            marginalGain[v_i]=prevD-d-dApprox;
            if (d+dApprox < CFGCC) {
              /*
              d+=computeExactDistance(E,reverse,dst,&dstApprox,Pinv);
              */
              CFGCC = d;
              tmpdst=dst;
              tmpApprox=dstApprox;
              s = v;
              bestMarginalGain=marginalGain[v_i];
            }
          }
        }
        S[i]=s;
        V[s]=false;
        Sdst=tmpdst;
        SdstApprox=tmpApprox;
        prevD=CFGCC;
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
      Lamg<CSRMatrix> lamg;
      CSRMatrix matrix = CSRMatrix::laplacianMatrix(G);
      lamg.setupConnected(matrix);
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
            for (count l = 0; l <vecOfSamples.size(); l++) {
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
    void CurrentFlowGroupCloseness::mergePeripheralNodes(){
      int pos;
      node c,s,w;
      double weight;
      std::vector<node> vecOfNodes,vecOfSupernodes;
      vecOfSupernodes.resize(0);
      vecOfNodes=G.nodes();
      for(count i=0;i<vecOfNodes.size();i++){
        c=vecOfNodes[i];
        if(G.degree(c)==1){
          s=G.randomNeighbor(c);
          auto it=std::find(vecOfSupernodes.begin(), vecOfSupernodes.end(), s);
          if(it== vecOfSupernodes.end()){
            vecOfSupernodes.push_back(s);
            vecOfPeripheralNodes.push_back(c);
          }
          else{
            pos = std::distance(vecOfSupernodes.begin(), it );
            w= vecOfPeripheralNodes[pos];
            weight=1./G.weight(s,w);
            weight+=1./G.weight(c,s);
            G.setWeight(s,w,1./weight);
            G.removeNode(c);
          }//end else
        }
      }
    }
    /***************************************************************************/
    void CurrentFlowGroupCloseness::coarseGraph(std::vector<node> vecOfChosenNodes,count degree){
      node c;
      double weight;
      std::vector<node> neighbors;
      std::vector<double> weights;

      weights.resize(degree);
      /*Path Coarsenening*/
      if(degree==2){
        for(count i=0;i<vecOfChosenNodes.size();i++){
          c=vecOfChosenNodes[i];
          neighbors=G.neighbors(c);
          weight=0.;
          for(count j=0;j<neighbors.size();j++){
            weights[j]=G.weight(c,neighbors[j]);
            weight+=weights[j];
          }
          /*Path Merge*/
          if(G.weight(neighbors[0],neighbors[1])!=0){
            weight=1./weight;
            weight+=1./G.weight(neighbors[0],neighbors[1]);
          }
          G.setWeight(neighbors[0],neighbors[1],1./weight);

          coarsedNodes.push_back(c);
          coarsedNeighbors.push_back(neighbors);
          coarsedWeights.push_back(weights);
          G.removeNode(c);
        }
      }
      else{
        for(count i=0;i<vecOfChosenNodes.size();i++){
          c=vecOfChosenNodes[i];
          computeStarCliqueWeights(c);
          neighbors=G.neighbors(c);
          for(count j=0;j<neighbors.size();j++){
            weights[j]=G.weight(c,neighbors[j]);
          }
          coarsedNodes.push_back(c);
          coarsedNeighbors.push_back(neighbors);
          coarsedWeights.push_back(weights);
          G.removeNode(c);
        }
      }
    }

    void CurrentFlowGroupCloseness::computeStarCliqueWeights(node c){
      node v,w,x;
      double Sum,weight;
      std::vector<node> vecOfNeighbors;

      vecOfNeighbors=G.neighbors(c);
      for(count i=0;i<vecOfNeighbors.size();i++){
        v=vecOfNeighbors[i];
        weight=1.;
        for(count j=0;j<vecOfNeighbors.size();j++){
          w=vecOfNeighbors[j];
          if(w!=v){
            weight*=G.weight(c,w);
          }
        }
        Sum+=weight;
      }
      for(count i=0;i<vecOfNeighbors.size();i++){
        v=vecOfNeighbors[i];
        for(count j=i+1;j<vecOfNeighbors.size();j++){
          w=vecOfNeighbors[j];
          weight=1.;
          for(count k=0;k<vecOfNeighbors.size();k++){
            x=vecOfNeighbors[k];
            if((x!=v)&&(x!=w)){
              weight*=G.weight(c,x);
            }//if x
          }//for x
          weight=Sum/weight;
          if(G.weight(v,w)!=0.){
            weight=1./weight;
            weight+=1./G.weight(v,w);
            weight=1./weight;
          }
          G.setWeight(v,w,weight);
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
      L= L+J;
      L=arma::inv(L);
      L= L-J;
      return L;
    }
    double CurrentFlowGroupCloseness::computeApproxDistances(std::vector<bool> E,std::vector<node> reverse, std::vector<double> dst, std::vector<double> *dstApprox){
      count w_j,l;
      node c,w;
      double d,distance;
      d=0.;
      std::vector<node> neighbors;
      std::vector<double> weights, tmp;

      tmp=*dstApprox;
      for(count i=0;i<coarsedNodes.size();i++){
        l=coarsedNodes.size()-i-1;
        neighbors=coarsedNeighbors[l];
        weights=coarsedWeights[l];
        distance=0.;
        for(count j=0;j<neighbors.size();j++){
          w=neighbors[j];
          w_j=reverse[w];
          if(E[w])
            distance=std::max(distance,dst[w_j]-weights[j]);
          else{
            distance=std::max(distance,tmp[w_j]-weights[j]);
          }
          tmp[l]=std::min(distance,tmp[l]);
        }
        d+=distance;
      }
      *dstApprox=tmp;
      return d;
    }
    double CurrentFlowGroupCloseness::computeExactDistance(std::vector<bool> E,std::vector<node> reverse,std::vector<double> dst, std::vector<double> *dstApprox, arma::Mat<double> Pinv){
      std::cout<<"Start EXACT"<<"\n";
      count l;
      double distance,distanceTMP,sum,weight,d;
      std::vector<count> neighbors_i;
      std::vector<double> cDistances,cDistancesTMP,xDistances,xDistancesTMP,weights,cdst,sdst,tmp;
      std::vector<node> neighbors;
      std::vector<std::vector<double>> coarsedDistances;
      arma::Mat<double> nDistances,nDistancesTMP;

      tmp=*dstApprox;
      coarsedDistances.resize(coarsedNodes.size());
      std::cout<<"Start CoarsedNodeLoop"<<"\n";
      d=0.;
      for(count i=0;i<coarsedNodes.size();i++){
        std::cout<<"Loop: "<< i<<"\n";
        l=coarsedNodes.size()-i-1;
        distance=0.;
        weights=coarsedWeights[l];
        neighbors=coarsedNeighbors[l];
        //Indices are computed central here
        neighbors_i.resize(neighbors.size());
        xDistances.resize(neighbors.size());
        std::cout<<"calc :"<< coarsedNodes[l] << "\n";
        std::cout<<"with neighbors: \n";
        for(count j=0;j<neighbors.size();j++){
          std::cout<< neighbors[j] << " , ";
          neighbors_i[j]=reverse[neighbors[j]];
          if(E[neighbors[j]])
            xDistances[j]=dst[neighbors_i[j]];
          else
            xDistances[j]=tmp[neighbors_i[j]];
        }
        std::cout<<"\n";
        std::cout<<"nDistances"<<"\n";
        nDistances=computeDistanceMatrix(neighbors,neighbors_i,E,Pinv);
        std::cout<<"cDistances"<<"\n";
        cDistances.resize(neighbors.size(),weights[0]);
        for(count j=1;j<neighbors.size();j++){
          cDistances[j]+=nDistances(0,j);
        }
        //Compute Sum
        sum=0;
        for(count a=0;a<neighbors.size();a++){
          weight=1.;
          for(count b=0;b<neighbors.size();b++){
            if(a!=b){
              weight*=weights[b];
            }
          }
          sum+=weight;
        }
        std::cout<<"FirstJoin"<<"\n";
        //FirstJoin
        distance=weights[0]+xDistances[0];
        std::cout<<"EdgeFire"<<"\n";
        //EdgeFire-Loop
        for(count j=1;j<weights.size();j++){
            distanceTMP=distance;
            xDistancesTMP=xDistances;
            cDistancesTMP=cDistances;
            nDistancesTMP=nDistances;
            distance-=residual(xDistancesTMP[j],distanceTMP,cDistancesTMP[j],0,cDistancesTMP[j],weights[j],1.);
          //Update xDistances after EdgeFire
            for(count a=0;a<neighbors.size();a++){
              xDistances[a]-=residual(xDistancesTMP[j],distanceTMP,nDistancesTMP(j,a),cDistancesTMP[a],cDistancesTMP[j],weights[j],1.);
              cDistances[a]-=residual(cDistancesTMP[j],0,nDistancesTMP(j,a),cDistancesTMP[a],cDistancesTMP[j],weights[j],1.);
            }
            //Update nDistances after EdgeFire
            for(count a=0;a<neighbors.size();a++){
              for(count b=a+1;b<neighbors.size();b++){
                nDistances(a,b)-=residual(nDistancesTMP(a,j),cDistances[a],nDistances(b,j),cDistances[b],cDistances[j],weights[j],1.);
                nDistances(b,a)=nDistances(a,b);
              }
            }
        }//end EDGE FIRE
        std::cout<<"Non-Bridge-Delete"<<"\n";
        //NonBrigdeDelete
      for(count a=0;a<neighbors.size();a++){
        for(count b=a+1;b<neighbors.size();b++){
          weight=1.;
          for(count j=0;j<weights.size();j++){
            if(j!=a&&j!=b){
              weight*=weights[j];
            }
          }
          distanceTMP=distance;
          xDistancesTMP=xDistances;
          cDistancesTMP=cDistances;
          nDistancesTMP=nDistances;
          distance+=residual(xDistancesTMP[b],xDistancesTMP[a],cDistancesTMP[b],cDistancesTMP[a],nDistancesTMP(a,b),weight,-1.);
        //Update xDistances after EdgeFire
          for(count y=0;y<neighbors.size();y++){
            xDistances[y]+=residual(xDistancesTMP[b],xDistancesTMP[a],nDistancesTMP(b,y),nDistancesTMP(a,y),nDistancesTMP(a,b),weight,-1.);
            cDistances[y]+=residual(cDistancesTMP[b],cDistancesTMP[a],nDistancesTMP(b,y),nDistancesTMP(a,y),nDistancesTMP(a,b),weight,-1);
          }
          //Update nDistances after EdgeFire
          for(count x=0;x<neighbors.size();x++){
            for(count y=0;y<neighbors.size();y++){
              nDistances(x,y)+=residual(nDistancesTMP(x,b),nDistancesTMP(x,a),nDistancesTMP(b,y),cDistances[b],nDistancesTMP(a,b),weight,-1.);
              nDistances(y,x)=nDistances(x,y);
            }
          }
        }
      }//End Non Bridge Delete
      std::cout<<"NonBrigdeDelete finished"<<"\n";
      coarsedDistances[l]=cDistances;
      tmp[l]=std::min(tmp[l],distance);
      d+=tmp[l];
      std::cout<<"save"<<"\n";
    }//End for(i)
      *dstApprox=tmp;
      std::cout<<"END EXACT"<<"\n";
      return d;
    }
    arma::Mat<double> CurrentFlowGroupCloseness::computeDistanceMatrix(std::vector<node> neighbors,std::vector<count> neighbors_i,std::vector<bool> E,arma::Mat<double> Pinv){
        arma::Mat<double> nDistances(neighbors.size(),neighbors.size());
        nDistances.zeros();
        for(count i=0;i<neighbors.size();i++){
          for(count j=i+1;j<neighbors.size();j++){
            if(E[neighbors[i]]&&E[neighbors[j]]){
              std::cout<<"case 1\n";
              nDistances(i,j)=Pinv(neighbors_i[i],neighbors_i[i])+Pinv(neighbors_i[j],neighbors_i[j])-2*Pinv(neighbors_i[i],neighbors_i[j]);
              nDistances(j,i)=nDistances(i,j);
            }
            else if(E[neighbors[i]]&&!E[neighbors[j]]){
              std::cout<<"case 2\n";
              nDistances(i,j)=distanceFinder(neighbors[i],neighbors_i[j]);
              nDistances(j,i)=nDistances(i,j);
            }
            else if(!E[neighbors[i]]&&E[neighbors[j]]){
              std::cout<<"case 3\n";
              std::cout<<"search for "<< neighbors[j] << " in " << neighbors[i] << "\n";
              nDistances(i,j)=distanceFinder(neighbors[j],neighbors_i[i]);
              nDistances(j,i)=nDistances(i,j);
            }
            else{
              std::cout<<"case 4\n";
              nDistances(i,j)=distanceFinder(coarsedNodes[std::max(neighbors_i[i],neighbors_i[j])],std::min(neighbors_i[i],neighbors_i[j]));
              nDistances(j,i)=nDistances(i,j);
            }
          }
        }
        std::cout<<"nDistances finished\n";
        return nDistances;
    }
    double CurrentFlowGroupCloseness::residual(double dstxj, double dstxi, double dstjy,double dstiy,double dstij,double weight,double factor){
      double variation;
      variation=dstxj-dstxi-dstjy+dstiy;
      variation*=variation;
      variation/=4*(weight+factor*dstij);
      return variation;
    }
    double CurrentFlowGroupCloseness::distanceFinder(node y, count x_i){
      std::cout<<"distance Finder\n";
      auto it=std::find(coarsedNeighbors[x_i].begin(), coarsedNeighbors[x_i].end(), y);
      int pos = std::distance(coarsedNeighbors[x_i].begin(), it );
      return coarsedDistances[x_i][pos];
      std::cout<<"distance Finder finished\n";
    }
} /* namespace NetworKit*/
