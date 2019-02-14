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

   CurrentFlowGroupCloseness::CurrentFlowGroupCloseness(Graph& G,const count k,const double epsilon,const count t): G(G),k(k),epsilon(epsilon),t(t){
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
    std::vector<node> CurrentFlowGroupCloseness::getNodesOfGroup(){
      return S;
    }
    /**
     * Returns maximum current flow group closeness
     */
    double CurrentFlowGroupCloseness::getCFGCC() {
      return CFGCC;
    }
    void CurrentFlowGroupCloseness::run(){
      if(t == 1){
        runDeterministicGreedy();
      }
      if(t == 2){
        runDeterministicGreedyLS();
      }
      if(t == 3){
        runTopK();
      }
      if(t == 4){
        runApproxJLT();
      }
      if(t == 5){
        runOptimum1();
      }
      if(t == 6 || t == 7){
        runRandom();
      }
      if(t == 8){
        runGCFCC2();
      }
    }

    void CurrentFlowGroupCloseness::runDeterministicGreedy(){
      count s, n;
      double d, max;
      std::vector<node> nodes, nbr;
      std::vector<double> gain;

      nodes=G.nodes();
      n=G.numberOfNodes();

      arma::Mat<double> L(G.numberOfNodes(),G.numberOfNodes());
      L.zeros();
      for(count i = 0; i < n; i++){
        nbr=G.neighbors(i);
        for(count j = 0; j < nbr.size(); j++){
          L(i,nbr[j]) = -1./G.weight(i,nbr[j]);
          L(i,i)+= 1./G.weight(i,nbr[j]);
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
      S[0]=s;
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

    void CurrentFlowGroupCloseness::runDeterministicGreedyLS(){
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
      std::cout << "2:"<< d << '\n';
      std::cout << "2:"<< (double)(n)/d << '\n';
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
        std::cout << "gain: " << max << '\n';
        S[l]=nodes[s];
        d-=max;
      }
     CFGCC = (double)(n)/d;
    }

    void CurrentFlowGroupCloseness::runApproxJLT(){
      count n=G.numberOfNodes();
      double max;
      std::vector<bool> av;
      std::vector<node> nodes = G.nodes();
      std::vector<count> re;
      std::vector<node> nbr;
      std::vector<double> gain;

      runApproxTopOne();

      count s = S[0];
      double d = CFGCC;
      std::cout << '\n';
      std::cout << "CFGCC: " << (double)(G.numberOfNodes())/ d << '\n';
      gain.resize(n,d);
      av.resize(G.upperNodeIdBound(),true);
      re.resize(G.upperNodeIdBound());
      epsilon = std::min(epsilon * d, 0.5);
      double ep = epsilon * epsilon;

      CSRMatrix L = CSRMatrix::laplacianMatrix(G);

      for(count l = 1; l < k; l++){

        L = shed(L, s);
        av[nodes[s]] = false;
        gain.erase(gain.begin()+s);
        nodes.erase(nodes.begin()+s);

        for(count i = 0; i < nodes.size(); i++){
          re[nodes[i]]=i;
        }

        n--;
        count q = std::min(ceil(4*log(n)/ep),(double)(n));
        double tol = (epsilon/(3*log(n)*log(n)))*sqrt((1-epsilon/3)/(1+epsilon/3));
        double randTab[2] = {-1/sqrt(q), 1/sqrt(q)};

        Lamg<CSRMatrix> lamg(tol);
        lamg.setupConnected(L);

        Vector x(n, 0.);
        std::vector<Vector> z1(q, Vector(n));
        std::vector<Vector> rhs1(q, Vector(n));
        std::vector<Vector> z2(q, Vector(n));
        std::vector<Vector> rhs2(q, Vector(n));
        std::vector<Vector> z3(q, Vector(n));
        std::vector<Vector> rhs3(q, Vector(n));

        for(count j = 0; j < l; j++){
          nbr = G.neighbors(S[j]);
          for(count w = 0; w < nbr.size(); w++){
            if(av[nbr[w]]){
              x[re[nbr[w]]] += 1.;
            }
          }
        }

        for(count i = 0; i < q; i++){

          G.forNodes([&](node u) {
            if(av[u]){
              rhs1[i][re[u]] = randTab[Aux::Random::integer(1)];
            }
          });

          G.forEdges([&](node u, node v) {
            double r = randTab[Aux::Random::integer(1)];
            if(av[u] & av[v]){
              if (u < v) {
                rhs2[i][re[u]] += r;
                rhs2[i][re[v]] -= r;
              }
              else {
                rhs2[i][re[u]] -= r;
                rhs2[i][re[v]] += r;
              }
            }
          });

          G.forNodes([&](node u) {
            if(av[u]){
              rhs3[i][re[u]] = randTab[Aux::Random::integer(1)] * sqrt(x[re[u]]);
            }
          });

        }//end q-loop

        lamg.parallelSolve(rhs1, z1);
        lamg.parallelSolve(rhs2, z2);
        lamg.parallelSolve(rhs3, z3);

        max=0;
        G.forNodes([&](node u) {
          double nume = 0.;
          double deno = 0.;
          if(av[u]){
            if(gain[re[u]]>max){
              for(count i = 0; i < q; i++){
                nume += pow(z1[i][re[u]],2);
                deno += pow(z2[i][re[u]],2) + pow(z3[i][re[u]],2);
              }
              gain[re[u]] = nume / deno;
              if (gain[re[u]] > max){
                max = gain[re[u]];
                s = re[u];
              }
            }
          }
        });
        std::cout << "gain: " << max << '\n';
        S[l] = nodes[s];
        av[nodes[s]] = false;
        d -= max;
      }
      CFGCC = (double)(G.numberOfNodes())/d;
    }

    void CurrentFlowGroupCloseness::runTopK(){
      std::vector<bool> s;
      count n=G.numberOfNodes();
      double d,max;
      std::vector<node> nbr;
      std::vector<double> scoreData;

      arma::Mat<double> L(G.numberOfNodes(),G.numberOfNodes());
      L.zeros();
      for(count i = 0; i < n; i++){
        nbr=G.neighbors(i);
        for(count j = 0; j < nbr.size(); j++){
          L(i,nbr[j]) = -1./G.weight(i,nbr[j]);
          L(i,i)+= 1./G.weight(i,nbr[j]);
        }
      }
      arma::Mat<double> J(L.n_rows,L.n_rows);
      J.fill(1./(double)(L.n_rows));
      arma::Mat<double> Inv= L+J;
      Inv = arma::inv_sympd(Inv);
      Inv = Inv-J;

      count index;
      bool search;
      s.resize(n,false);
      for(count l = 0; l < k; l++){
        index = 0;
        while(search) {
          if(s[index])
            search = false;
        }
        for (count i = 0; i < n; i++) {
          if(!s[i]){
            if(Inv(i,i)< Inv(index,index)){
              d = scoreData[i];
              index = i;
            }
          }
        }
        s[index]=true;
      }
      count j=0;
      for (count i = 0; i < n; i++) {
        if(s[i]){
          L.shed_row(s[i]);
          L.shed_col(s[i]);
          S[j]=i;
          j++;
        }
      }
     L=arma::inv_sympd(L);

     CFGCC = (double)(n)/arma::trace(L);
    }

    void CurrentFlowGroupCloseness::runApproxTopOne(){
      const count n = G.numberOfNodes();
      double ep = pow(epsilon/3.,2)/2 + pow(epsilon/3.,3)/3;
      const count q = std::min(ceil(4*log(n)/ep),(double)(n));
      const double tol = (epsilon/(3*n*n))*sqrt((1-epsilon/3)/(1+epsilon/3));
      double randTab[2] = {1/sqrt(q), -1/sqrt(q)};
      std::vector<double> scoreData;

      scoreData.resize(n, 0.);
      std::vector<Vector> rhs(q, Vector(n));
      std::vector<Vector> z(q, Vector(n));

      Lamg<CSRMatrix> lamg(tol);
      CSRMatrix L = CSRMatrix::laplacianMatrix(G);
      lamg.setupConnected(L);

      for (index i = 0; i < q; i++) {
        G.forEdges([&](node u, node v) {
          double r = randTab[Aux::Random::integer(1)];
          if (u < v) {
            rhs[i][u] += r;
            rhs[i][v] -= r;
          }
          else {
            rhs[i][u] -= r;
            rhs[i][v] += r;
          }
        });
      }

      lamg.parallelSolve(rhs, z);

      for (index i = 0; i < q; i++){
        G.forNodes([&](node u){
          double z2 = z[i][u];
          G.forNodes([&](node v){
            double z3 = fabs(z2 - z[i][v]);
            scoreData[u] += z3*z3;
          });
        });
      }

      double max=0;
      node s = G.randomNode();
      G.forNodes([&](node u){
        if(scoreData[u]<scoreData[s]){
          s=u;
        }
      });

      S[0] = s;
      CFGCC = scoreData[s];

    }

    void CurrentFlowGroupCloseness::runOptimum1(){
      double COPT;
      std::vector<bool> s;
      std::vector<node> SOPT;

      COPT=0;
      s.resize(G.numberOfNodes(),false);
      for (count i = 0; i < k; i++) {
        s[i]=true;
      }
      do{
        closeness1(s);
        if(CFGCC>COPT){
          COPT=CFGCC;
          SOPT=S;
        }
      }while(std::next_permutation(s.begin(), s.end()));

      CFGCC=COPT;
      S=SOPT;
    }


    void CurrentFlowGroupCloseness::runRandom(){
      double p;
      std::vector<bool> s;

      s.resize(G.numberOfNodes(),false);

      for (count i = 0; i < k; i++) {
        p = Aux::Random::integer(G.numberOfNodes());
        while(s[p] == true){
          p = Aux::Random::integer(G.numberOfNodes());
        }
        s[p] = true;
      }

      if(t == 6){closeness1(s);}
      if(t == 7){closeness2(s);}
    }

    void CurrentFlowGroupCloseness::closeness1(std::vector<bool> s){
      count n;
      double d;
      std::vector<node> nbr;

      n=G.numberOfNodes();

      arma::Mat<double> L(G.numberOfNodes(),G.numberOfNodes());
      L.zeros();
      for(count i = 0; i < n; i++){
        nbr=G.neighbors(i);
        for(count j = 0; j < nbr.size(); j++){
          L(i,nbr[j]) = -1./G.weight(i,nbr[j]);
          L(i,i)+= 1./G.weight(i,nbr[j]);
        }
      }
      count j=0;
      for (count i = 0; i < n; i++) {
        if(s[i]){
          L.shed_row(i-j);
          L.shed_col(i-j);
          S[j]=i;
          j++;
        }
      }
     L=arma::inv_sympd(L);

     CFGCC = (double)(n)/arma::trace(L);
    }

    void CurrentFlowGroupCloseness::closeness2(std::vector<bool> s){
      count n;
      double d, m;
      std::vector<node> nbr;
      n=G.numberOfNodes();

      arma::Mat<double> L(G.numberOfNodes(),G.numberOfNodes());
      L.zeros();
      for(count i = 0; i < n; i++){
        nbr=G.neighbors(i);
        for(count j = 0; j < nbr.size(); j++){
          L(i,nbr[j]) = -1./G.weight(i,nbr[j]);
          L(i,i)+= 1./G.weight(i,nbr[j]);
        }
      }
      arma::Mat<double> J(L.n_rows,L.n_rows);
      J.fill(1./(double)(L.n_rows));
      L = L+J;
      L = arma::inv_sympd(L);
      L = L - J;

      d = 0.;
      for(count i = 0; i < n; i++){
        m = 2 * L.max();
        for (count j = 0; j < n; j++) {
          if(s[j]){
            if((L(i,i)+L(j,j)-2*L(i,j))<m){
              m = L(i,i)+L(j,j)-2*L(i,j);
            }
          }
        }
        d += m;
      }
      CFGCC = (double)(n)/d;
    }

    void CurrentFlowGroupCloseness::runGCFCC2(){
      count n;
      node tmps;
      double m, r, r2, tmpr;
      std::vector<node> nbr;
      std::vector<double> d, d2, tmpd, g;
      std::vector<bool> s;
      n=G.numberOfNodes();

      arma::Mat<double> L(G.numberOfNodes(),G.numberOfNodes());
      L.zeros();
      for(count i = 0; i < n; i++){
        nbr=G.neighbors(i);
        for(count j = 0; j < nbr.size(); j++){
          L(i,nbr[j]) = -1./G.weight(i,nbr[j]);
          L(i,i)+= 1./G.weight(i,nbr[j]);
        }
      }
      arma::Mat<double> J(L.n_rows,L.n_rows);
      J.fill(1./(double)(L.n_rows));
      L = L+J;
      L = arma::inv_sympd(L);
      L = L - J;
      d.resize(n,2 * L.max());
      g.resize(n,n * 2 * L.max());
      s.resize(n,true);
      r = n * 2 * L.max();
      for(count l = 0; l < k; l++){
        m = 0;
        r2 = r;
        for(count i = 0; i < n; i++){
          if(s[i] && (g[i] >= m)){
            tmpd = d;
            tmpr = 0.;
            for (count j = 0; j < n; j++) {
              if((L(i,i)+L(j,j)-2*L(i,j))<d[j]){
                tmpd[j] = L(i,i)+L(j,j)-2*L(i,j);
              }
              tmpr += tmpd[j];
            }
            g[i] = r2 - tmpr;
            if(tmpr < r){
              tmps = i;
              r = tmpr;
              m = g[i];
              d2 = tmpd;
            }
          }
        }
        d = d2;
        s[tmps] = false;
      }
      count j = 0;
      std::vector<node> nodes = G.nodes();
      for(count i = 0; i < n; i++){
        if(!s[i]){
          S[j] = nodes[i];
          j++;
        }
      }
      CFGCC = (double)(n)/r;
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
