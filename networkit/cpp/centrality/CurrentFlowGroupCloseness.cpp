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
#include "../numerics/ConjugateGradient.h"
#include "../auxiliary/Log.h"
#include "../numerics/Preconditioner/DiagonalPreconditioner.h"
#include <chrono>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include "../components/ConnectedComponents.h"
#include <armadillo>

namespace NetworKit {

   CurrentFlowGroupCloseness::CurrentFlowGroupCloseness(Graph& G,const count k,const double epsilon): G(G),k(k),epsilon(epsilon){
     if (G.isDirected()) throw std::runtime_error("Graph is directed!");
     ConnectedComponents cc(G);
     cc.run();
     if (cc.getPartition().numberOfSubsets() > 1) throw std::runtime_error("Graph has more then one component!");
     if(k>=G.numberOfNodes()) throw std::runtime_error("Size of Group greater then number of nodes!");

     S.resize(k);
     CFGCC=0.;
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
    void CurrentFlowGroupCloseness::runExact(){
      count s, n;
      double d, max;
      std::vector<node> nodes, nbr;
      std::vector<count> re;
      std::vector<double> gain;

      nodes = G.nodes();
      n=G.numberOfNodes();

      re.resize(G.upperNodeIdBound());
      for(count i = 0; i < nodes.size(); i++){
        re[nodes[i]]=i;
      }
      arma::Mat<double> L(G.numberOfNodes(),G.numberOfNodes());
      L.zeros();
      for(count i = 0; i < nodes.size(); i++){
        nbr=G.neighbors(nodes[i]);
        for(count j = 0; j < nbr.size(); j++){
          L(i,re[nbr[j]]) = -1./G.weight(nodes[i],nbr[j]);
          L(i,i)+= 1./G.weight(nodes[i],nbr[j]);
        }
      }
      arma::Mat<double> J(L.n_rows,L.n_rows);
      J.fill(1./(double)(L.n_rows));
      arma::Mat<double> Inv= L+J;
      Inv = arma::inv_sympd(Inv);
      Inv = Inv-J;

      s = 0;
      d = Inv(0,0);
      for (count i = 1; i < n; i++) {
        if(Inv(i,i)< d){
          d = Inv(i,i);
          s = i;
        }
      }
      S[0]=nodes[s];
      L.shed_row(s);
      L.shed_col(s);
      nodes.erase(nodes.begin()+s);
      Inv=arma::inv(L);

      gain.resize(nodes.size(),arma::trace(Inv));
      for(count l = 1; l < k; l++){
        max=0.;
        for (count i = 0; i < nodes.size(); i++) {
          if(gain[i]>max){
            gain[i]=(Inv.row(i)*Inv.col(i)).eval()(0,0)/Inv(i,i);
            if (gain[i] > max){
              max=gain[i];
              s=i;
            }
          }
        }
        S[l]=nodes[s];
        J=Inv;
        for(count i = 0; i < nodes.size(); i++){
          for(count j = i; j < nodes.size(); j++){
            Inv(i,j)-=J(s,i)*J(s,j)/J(s,s);
            Inv(j,i)=Inv(i,j);
          }
        }
        Inv.shed_row(s);
        Inv.shed_col(s);
        nodes.erase(nodes.begin()+s);
        gain.erase(gain.begin()+s);
      }
     CFGCC = (double)(n)/arma::trace(Inv);
    }

    void CurrentFlowGroupCloseness::runExactLS(){
      node s;
      count n;
      double d,max;
      std::vector<node> nodes;
      std::vector<count> shed_indices;
      std::vector<double> gain;

      n=G.numberOfNodes();
      shed_indices.resize(1,0);
      ConjugateGradient<CSRMatrix,DiagonalPreconditioner> linearSolver;
      CSRMatrix matrix = CSRMatrix::laplacianMatrix(G);
      matrix = shed( matrix, 0);
      linearSolver.setupConnected(matrix);

      nodes = G.nodes();

      gain.resize(nodes.size(),0.);
      Vector result(nodes.size()-1);
      Vector rhs(nodes.size()-1, 0.);
      Vector zero(nodes.size()-1, 0.);
      for (count i = 0; i < n-1; i++) {
        result=zero;
        rhs[i]=1.;
        linearSolver.solve(rhs, result);
        gain[0]+=result[i];
        for(count j = 0; j < n-1; j++){
          gain[i+1]+=result[i]-2.*result[j];
          gain[j+1]+=result[i];
        }
        rhs[i]=0.;
      }
      s = 0;
      d = gain[0];
      for (count i = 1; i < n; i++) {
        if(gain[i]< d){
          d = gain[i];
          s = i;
        }
      }

      S[0]=nodes[s];
      matrix=CSRMatrix::laplacianMatrix(G);

      for(count l = 1; l < k; l++){
        matrix = shed(matrix, s);
        gain.erase(gain.begin()+s);
        nodes.erase(nodes.begin()+s);
        linearSolver.setupConnected(matrix);
        Vector result(nodes.size(), 0.);
        Vector rhs(nodes.size(), 0.);
        Vector zero(nodes.size(), 0.);
        max=0.;
        for (count i = 0; i < nodes.size(); i++) {
          if(gain[i]>max){
            result=zero;
            rhs[i]=1.;
            linearSolver.solve(rhs, result);
            gain[i]=result*result;
            gain[i]/=result[i];
            rhs[i]=0.;
            if (gain[i] > max){
              max=gain[i];
              s=i;
            }
          }
        }
        S[l]=nodes[s];
        d-=max;
      }
     CFGCC = (double)(n)/d;
    }

    void CurrentFlowGroupCloseness::runApproxJLT(){
      node s;
      count n;
      double d,max;
      std::vector<node> nodes;
      std::vector<double> gain;

      n=G.numberOfNodes();
      shed_indices.resize(1,0);
      ConjugateGradient<CSRMatrix,DiagonalPreconditioner> linearSolver;
      CSRMatrix matrix = CSRMatrix::laplacianMatrix(G);
      matrix = shed( matrix, 0);
      linearSolver.setupConnected(matrix);

      nodes = G.nodes();

      gain.resize(nodes.size(),0.);
      Vector result(nodes.size()-1);
      Vector rhs(nodes.size()-1, 0.);
      Vector zero(nodes.size()-1, 0.);
      for (count i = 0; i < n-1; i++) {
        result=zero;
        rhs[i]=1.;
        linearSolver.solve(rhs, result);
        gain[0]+=result[i];
        for(count j = 0; j < n-1; j++){
          gain[i+1]+=result[i]-2.*result[j];
          gain[j+1]+=result[i];
        }
        rhs[i]=0.;
      }
      s = 0;
      d = gain[0];
      for (count i = 1; i < n; i++) {
        if(gain[i]< d){
          d = gain[i];
          s = i;
        }
      }

      S[0]=nodes[s];
      matrix=CSRMatrix::laplacianMatrix(G);

      for(count l = 1; l < k; l++){
        matrix = shed(matrix, s);
        gain.erase(gain.begin()+s);
        nodes.erase(nodes.begin()+s);
        linearSolver.setupConnected(matrix);
        Vector result(nodes.size(), 0.);
        Vector rhs(nodes.size(), 0.);
        Vector zero(nodes.size(), 0.);
        max=0.;
        for (count i = 0; i < nodes.size(); i++) {
          if(gain[i]>max){
            result=zero;
            rhs[i]=1.;
            linearSolver.solve(rhs, result);
            gain[i]=result*result;
            gain[i]/=result[i];
            rhs[i]=0.;
            if (gain[i] > max){
              max=gain[i];
              s=i;
            }
          }
        }
        S[l]=nodes[s];
        d-=max;
      }
     CFGCC = (double)(n)/d;
    }



    CSRMatrix CurrentFlowGroupCloseness::shed(CSRMatrix matrix, count s){
      count j=0;
      CSRMatrix left=CSRMatrix(matrix.numberOfRows()-1,matrix.numberOfColumns());
      CSRMatrix right=CSRMatrix(matrix.numberOfRows(),matrix.numberOfColumns()-1);
      CSRMatrix sub=CSRMatrix(matrix.numberOfRows()-1,matrix.numberOfColumns()-1);
      for(count i = 0; i < matrix.numberOfRows();i++){
        if(i!=s){
          left.setValue(j,i,1.);
          right.setValue(i,j,1.);
          j++;
        }
      }
      sub = left * matrix;
      sub = sub * right;
      return sub;
    }
} /* namespace NetworKit*/
