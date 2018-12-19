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
     E.resize(limit,true);

     end = std::chrono::high_resolution_clock::now();
     diff = end-start;
     std::cout << "Constructor finished in " << diff.count() << "(s)" << "\n";
     std::cout << "Number of Nodes=" << G.numberOfNodes() << "\n";
   }

   void CurrentFlowGroupCloseness::run() {
     bool coarse;
     count degree;
     auto start = std::chrono::high_resolution_clock::now();
     auto end = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double> diff;
     std::vector<node> chosenNodes;

     degree=minDegree();
     coarse=true;
     while((degree<CB)&&(coarse)){
       chosenNodes=coarsingNodes(degree, false);
       if(G.numberOfNodes()-chosenNodes.size()>k){
          graphCoarsening(chosenNodes);
          degree=minDegree();
        }
       else
         coarse=false;
     }

     end = std::chrono::high_resolution_clock::now();
     diff = end-start;
     std::cout << "Coarsening of G finished in: " << diff.count() << "(s)" << "\n";
     std::cout << "Nodes: " << G.numberOfNodes() << ", Coarsed: "<< coarsedNodes.size()<<"("<<(double)(coarsedNodes.size())/(double)(n)*100.<<"%)\n";
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
    /**
     * Returns group of size k with maximum closeness.
     */
    std::vector<node> CurrentFlowGroupCloseness::getNodesofGroup(){
      return S;
    }
    /**
     * Returns maximum current flow group closeness
     */
    double CurrentFlowGroupCloseness::getCFGCC() {
      return CFGCC;
    }
    /**
     *
     */
    void CurrentFlowGroupCloseness::greedy(){
      bool search;
      count a;
      node tmps,v,w;
      double sd,d,min,exact, tmpd, tmpGain, distance;
      std::vector<bool> V, E;
      std::vector<node> nodes;
      std::vector<double> sdst, dst, tmpdst, gain;
      arma::Mat<double> Pinv;

      tmpd=n*n*n;
      sd=tmpd;
      V.resize(G.upperNodeIdBound(),true);
      sdst.resize(limit,n*n);
      gain.resize(limit,tmpd);

      nodes = G.nodes();
      coarsedDistances.resize(coarsedNodes.size());

      Pinv=pinv();
      for(count l = 0; l < k; l++){
        tmpGain = 0.;
        for (count i = 0; i < nodes.size(); i++) {
          if(V[nodes[i]] && (tmpGain<gain[nodes[i]])){
            //reset distances
            dst = sdst;
            d = 0.;
            min = 0.;
            //Existing nodes
            for (count j = 0; j < nodes.size(); j++) {
              distance = Pinv(i,i) + Pinv(j,j) - 2 * Pinv(i,j);
              dst[nodes[j]]=std::min(dst[nodes[j]],distance);
              d += dst[nodes[j]];
            }
            //Min distances to coarsed nodes
            for(count j = 0; j < coarsedNodes.size(); j++){
              a = coarsedNodes.size()-j-1;
              distance = 0;
              for(count b = 0; b < coarsedNeighbors[a].size(); b++){
                distance = std::max(distance, dst[coarsedNeighbors[a][b]] - coarsedWeights[a][b]);
              }
              min += std::min(dst[coarsedNodes[a]],distance);
            }
            gain[nodes[i]] = sd - d - min;
            if (d+min< tmpd) {
              exact=0;
              search=true;
              a = coarsedNodes.size()-1;
              std::cout << "/* computing node */"<< nodes[i] << '\n';
              while(search){
                distance=exactDistance(a,dst,Pinv);
                dst[coarsedNodes[a]] = std::min(dst[coarsedNodes[a]], distance);
                exact += dst[coarsedNodes[a]];
                std::cout << "/* a */"<< a << '\n';
                if((d+exact>=tmpd)||(a==0))
                  search=false;
                a--;
              }
              gain[nodes[i]] = sd - d - exact;
              if(d+exact<tmpd){
                tmpd = d + exact;
                tmpdst = dst;
                tmpGain = gain[nodes[i]];
                tmps = nodes[i];
              }
            }
          }
        }
        /*
        for(count i=0; i<coarsedNodes.size();i++){
          for(count j=0;j<nodes();j++)
        }
        */
        S[l]=tmps;
        V[tmps]=false;
        sdst=tmpdst;
        sd=tmpd;
        std::cout << "S["<<l<<"]: " << S[l] << " sd: "<< sd << '\n';
      }
     CFGCC = (double)(n)/sd;
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
        vecOfSamples.push_back(v);
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
    /*
    *
    *returns min degree of graph G
    */
    count CurrentFlowGroupCloseness::minDegree(){
      bool search=true;
      count min=CB;
      count i=0;
      std::vector<node> nodes;

      nodes=G.nodes();

      while((i<nodes.size())&&(search)){
        if(G.degree(nodes[i])<min){
          min=G.degree(nodes[i]);
          if(min==1){
            search=false;
          }
        }
        i++;
      }
      return min;
    }
    /*
    *
    *@param degree
    *@param Random
    *
    *creates a vector of nodes with input degree, in the list no pair of
    *nodes share a same neighbor
    */
    std::vector<node> CurrentFlowGroupCloseness::coarsingNodes(count degree, bool Random){
        bool search;
        count l;
        std::vector<bool> A;
        std::vector<node> nodes, chosenNodes, neighbors, candidates;

        chosenNodes.resize(0);
        candidates.resize(0);
        nodes=G.nodes();
        A.resize(G.upperNodeIdBound(),true);

        for(count i=0;i<nodes.size();i++){
          if(G.degree(nodes[i])==degree){
            candidates.push_back(nodes[i]);
          }
        }
        if(Random){
          std::random_shuffle (candidates.begin(), candidates.end());
        }

        for(count i=0;i<candidates.size();i++){
          if(A[candidates[i]]){
            neighbors=G.neighbors(candidates[i]);
            search=true;
            l=0;
            if((search)&&(l<neighbors.size())){
              if(!A[neighbors[l]]){
                search=false;
              }
              else{
                l++;
              }
            }
            if(search){
              for(count j=0;j<neighbors.size();j++){
                A[neighbors[j]]=false;
              }
              A[candidates[i]]=true;
              chosenNodes.push_back(candidates[i]);
            }
          }
        }
        return chosenNodes;
      }
    /**
     *
     * @param nodes, vector of nodes that will be coarsed in this round
     *
     *Computes via adjacentCliqueWeights corresponding weights that are necessary
     *to remain distances and removes all nodes in vector nodes and
     *stores: coarsed nodes with corresponding weights and neighbors
     */
    void CurrentFlowGroupCloseness::graphCoarsening(std::vector<node> nodes){
      std::vector<node> neighbors;
      std::vector<double> weights;

      weights.resize(G.degree(nodes[0]));

      for(count i=0;i<nodes.size();i++){
        neighbors=G.neighbors(nodes[i]);
        for(count j=0;j<neighbors.size();j++){
          weights[j]=G.weight(nodes[i],neighbors[j]);
        }
        adjacentCliqueWeights(neighbors, weights);
        coarsedNodes.push_back(nodes[i]);
        coarsedNeighbors.push_back(neighbors);
        coarsedWeights.push_back(weights);
        G.removeNode(nodes[i]);
        E[nodes[i]]=false;
      }
    }
    /**
     *
     * @param neighbors, vector of nodes that will joined as a clique
     * @param weights, will be added to the clique
     *
     *Computes the clique corresponding to c in function graphCoarsening in two steps
     *First: precomputes sum, Second: calculates for each edge pair the weight
     *and adds the edge, (if necessary it merges two edges)
     */
    void CurrentFlowGroupCloseness::adjacentCliqueWeights(std::vector<node> neighbors,std::vector<double> weights){
      double sum,weight;

      for(count i=0;i<neighbors.size();i++){
        weight=1.;
        for(count l=0;l<neighbors.size();l++){
          if(i!=l){
            weight*=weights[l];
          }
        }
        sum+=weight;
      }
      for(count i=0;i<neighbors.size();i++){
        for(count j=i+1;j<neighbors.size();j++){
          weight=1.;
          for(count l=0;l<neighbors.size();l++){
            if((l!=i)&&(l!=j)){
              weight*=weights[l];
            }//if x
          }//for x
          weight=sum/weight;
          if(G.weight(neighbors[i],neighbors[j])!=0.){
            weight=1./weight;
            weight+=1./G.weight(neighbors[i],neighbors[j]);
            weight=1./weight;
          }
          G.setWeight(neighbors[i],neighbors[j],weight);
        }// for w
      }// for v
    }
    /*
    *
    *Returns Pseudo Inverse of Laplacian.
    *
    *First:Creates a Laplacian for G
    *Second: Computes a pertubation and Invertes it.
    */
    arma::Mat<double> CurrentFlowGroupCloseness::pinv(){
      double factor, weight;
      std::vector<node> nodes, nbr;
      std::vector<count> re;
      arma::Mat<double> L(G.numberOfNodes(),G.numberOfNodes());
      L.zeros();

      nodes=G.nodes();

      re.resize(G.upperNodeIdBound());
      for(count i = 0; i<nodes.size(); i++){
        re[nodes[i]]=i;
      }

      for(count i = 0; i < nodes.size(); i++){
        nbr=G.neighbors(nodes[i]);
        for(count j = 0; j < nbr.size(); j++){
          weight = 1./G.weight(nodes[i],nbr[j]);
          L(i,re[nbr[j]]) = -weight;
          L(i,i)+= weight;
        }
      }

      arma::Mat<double> J(L.n_rows,L.n_rows);
      factor=1./(double)(L.n_rows);
      J.fill(factor);
      L= L+J;
      L=arma::inv_sympd(L);
      L= L-J;
      return L;
    }
    /*
    *
    */
    double CurrentFlowGroupCloseness::exactDistance(count c,std::vector<double> dst, arma::Mat<double> Pinv){
      double distance, distanceTMP, sum, weight, factor, variation;
      std::vector<double> cDistances,cDistancesTMP,xDistances,xDistancesTMP,weights;
      std::vector<node> nbr;
      arma::Mat<double> nDistances,nDistancesTMP;

      weights=coarsedWeights[c];
      nbr=coarsedNeighbors[c];
      xDistances.resize(nbr.size());
      cDistances.resize(nbr.size());

      //Compute xDistances
      for(count i=0;i<nbr.size();i++){
        xDistances[i]=dst[nbr[i]];
      }
      //Compute nDistances extern
      nDistances=distanceMatrix(nbr,Pinv);

      //First Join of c <-> Compute cDistances
      for(count i=0;i<nbr.size();i++){
        cDistances[i]+=nDistances(0,i)+weights[i];
      }

      //Distance after First Join
      distance=xDistances[0]+cDistances[0];
      //EdgeFire-Loop
      for(count i=1;i<weights.size();i++){
        //Store TMP Values
        distanceTMP=distance;
        xDistancesTMP=xDistances;
        cDistancesTMP=cDistances;
        nDistancesTMP=nDistances;
        //Factor
        factor=4*(weights[i]+cDistancesTMP[i]);
        /*
        * Residual(1.d(x,i)-2.d(x,j)-3.d(i,y)-4.d(j,y)-5.d(i,j),6.weight(i,j))
        * For distance x=x, y=c, i=i, j=c
        * => d(x,i)=xDistance[i], d(x,j)=d(x,c), d(i,y)=d(c,i), d(j,y)=0., d(i,j)=d(c,i)
        */
        distance-=residual(xDistancesTMP[i]-distanceTMP-cDistancesTMP[i],factor);
        /*
        * For xDistances x=x, y=y, i=i, j=c
        */
        variation=xDistancesTMP[i]-distanceTMP;
        for(count y=0;y<nbr.size();y++){
          xDistances[y]-=residual(variation-nDistancesTMP(i,y)+cDistancesTMP[y], factor);
        }
        /*
        * For cDistances x=c, y=y, i=i, j=c
        */
        variation=cDistancesTMP[i];
        for(count y=0;y<nbr.size();y++){
          cDistances[y]-=residual(variation-nDistancesTMP(i,y)+cDistancesTMP[y],factor);
        }
        /*
        * Residual(1.d(x,i)-2.d(x,j)-3.d(i,y)-4.d(j,y)-5.d(i,j),6.weight(i,j))
        * For nDistances x=a, y=b, i=i, j=c
        */
        for(count a=0;a<nbr.size();a++){
          variation=nDistancesTMP(a,i)-cDistances[a];
          for(count b=a+1;b<nbr.size();b++){
            nDistances(a,b)-=residual(variation-nDistances(i,b)+cDistances[b],factor);
            nDistances(b,a)=nDistances(a,b);
          }
        }
      }//end EDGE FIRE

      //NonBrigdeDelete
      //Compute sum
      sum=starWeight(c);
      //NonBridgeDelete Loop
      for(count a=0;a<nbr.size();a++){
        for(count b=a+1;b<nbr.size();b++){
          weight=1.;
          for(count j=0;j<weights.size();j++){
            if( (j!=a) && (j!=b) ){
              weight *= weights[j];
            }
          }
          weight=sum/weight;
          //Store TMP values
          distanceTMP=distance;
          xDistancesTMP=xDistances;
          cDistancesTMP=cDistances;
          nDistancesTMP=nDistances;
          //
          factor=4*(weight-nDistances(a,b));
          /*
          * Residual(1.d(x,i)-2.d(x,j)-3.d(i,y)-4.d(j,y)-5.d(i,j),6.weight(i,j))
          * For distance x=x, y=c, i=a, j=b
          * => d(x,i)=xDistance[a], d(x,j)=d(x,a), d(i,y)=d(c,a), d(b,y)=0., d(i,j)=d(a,b)
          */

          distance += residual(xDistancesTMP[a]-xDistancesTMP[b]-cDistancesTMP[a]+cDistancesTMP[b],factor);
          /*
          * Residual(1.d(x,i)-2.d(x,j)-3.d(i,y)-4.d(j,y)-5.d(i,j),6.weight(i,j))
          * For xDistances x=x, y=y, i=a, j=b
          */
          variation = xDistancesTMP[a] - xDistancesTMP[b];
          for(count y=0;y<nbr.size();y++){
            xDistances[y] += residual(variation-nDistancesTMP(a,y)+nDistancesTMP(b,y),factor);
          }
          /*
          * For cDistances c=x, y=y, i=a, j=b
          */
          variation = cDistancesTMP[a] - cDistancesTMP[b];
          for(count y=0;y<nbr.size();y++){
            cDistances[y] += residual(variation - nDistancesTMP(a,y) + nDistancesTMP(b,y),factor);
          }
          /*
          * Residual(1.d(x,i)-2.d(x,j)-3.d(i,y)-4.d(j,y)-5.d(i,j),6.weight(i,j))
          * For nDistances x=x, y=y, i=a, j=b
          */
          for(count x=0;x<nbr.size();x++){
            variation = nDistancesTMP(x,a) - nDistancesTMP(x,b);
            for(count y=x+1;y<nbr.size();y++){
              nDistances(x,y)+=residual(variation - nDistancesTMP(a,y) + nDistancesTMP(b,y),factor);
              nDistances(y,x)=nDistances(x,y);
            }
          }
        }
      }//End Non Bridge Delete
      //Store Global!!
      coarsedDistances[c]=cDistances;
      return distance;
    }
    /*
    *
    */
    arma::Mat<double> CurrentFlowGroupCloseness::distanceMatrix(std::vector<node> nbr,arma::Mat<double> Pinv){
        std::cout << "nDistances start" << '\n';
        std::vector<count> re;
        std::vector<node> nodes;
        arma::Mat<double> nDistances(nbr.size(),nbr.size());
        nDistances.zeros();

        nodes=G.nodes();
        re.resize(G.upperNodeIdBound());
        for(count i = 0; i<nodes.size(); i++){
          re[nodes[i]]=i;
        }
        for(count i=0; i<coarsedNodes.size(); i++){
          re[coarsedNodes[i]]=i;
        }

        for(count i=0;i<nbr.size();i++){
          for(count j=i+1;j<nbr.size();j++){
            if(E[nbr[i]]&&E[nbr[j]]){
              nDistances(i,j)=Pinv(re[nbr[i]],re[nbr[i]])+Pinv(re[nbr[j]],re[nbr[j]])-2*Pinv(re[nbr[i]],re[nbr[j]]);
              nDistances(j,i)=nDistances(i,j);
            }
            else if(E[nbr[i]]&&!E[nbr[j]]){
              nDistances(i,j)=distanceFinder(nbr[i],re[nbr[j]]);
              nDistances(j,i)=nDistances(i,j);
            }
            else if(!E[nbr[i]]&&E[nbr[j]]){
              nDistances(i,j)=distanceFinder(nbr[j],re[nbr[i]]);
              nDistances(j,i)=nDistances(i,j);
            }
            else{
              nDistances(i,j)=distanceFinder(coarsedNodes[std::max(re[nbr[i]],re[nbr[j]])],std::min(re[nbr[i]],re[nbr[j]]));
              nDistances(j,i)=nDistances(i,j);
            }
          }
        }
        std::cout << "nDistances finished" << '\n';
        return nDistances;
    }
    /*
    *
    */
    double CurrentFlowGroupCloseness::residual(double variation, double factor){
      variation*=variation;
      variation/=factor;
      return variation;
    }
    /*
    *
    */
    double CurrentFlowGroupCloseness::distanceFinder(node y, count x_i){
      auto it=std::find(coarsedNeighbors[x_i].begin(), coarsedNeighbors[x_i].end(), y);
      int pos = std::distance(coarsedNeighbors[x_i].begin(), it );
      return coarsedDistances[x_i][pos];
    }
    /*
    *Computes Sum as Clique-Weight-Factor
    */
    double CurrentFlowGroupCloseness::starWeight(count c){
      double sum=0;
      double weight;

      for(count i=0;i<coarsedNeighbors[c].size();i++){
        weight=1.;
        for(count j=0;j<coarsedNeighbors[c].size();j++){
          if(i!=j){
            weight*=coarsedWeights[c][j];
          }
        }
        sum+=weight;
      }
      return sum;
    }
} /* namespace NetworKit*/
