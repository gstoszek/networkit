/*
*  EffectiveResistanceDistance.cpp
*
*      Author: gstoszek
*/

#include "EffectiveResistanceDistance.h"
#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include "../auxiliary/Log.h"
#include <chrono>
#include <stdlib.h>
#include <cmath>

namespace NetworKit {


  EffectiveResistanceDistance::EffectiveResistanceDistance(const count n){
    M.set_size(n,n);
    M.zeros();
  };

  void EffectiveResistanceDistance::computeFromLaplacian(std::vector<node> vecOfNodes,arma::Mat<double> Laplacian){
    node v;
    node w;
    count n;
    double factor;
    n=Laplacian.n_rows;
    arma::Mat<double> J(Laplacian.n_rows,Laplacian.n_rows);
    factor=1./n;
    J.fill(factor);
    Laplacian= Laplacian+J;
    Laplacian=arma::inv_sympd(Laplacian);
    Laplacian= Laplacian-J;
    /*calculate initial M Matrix*/
    for(count i=0;i<vecOfNodes.size();i++){
      v=vecOfNodes[i];
      for(count j=i+1;j<vecOfNodes.size();j++){
        w=vecOfNodes[j];
        M(v,w)=Laplacian(i,i)+Laplacian(j,j)-2.*Laplacian(i,j);
        M(w,v)=M(v,w);
      }
    }
  }
  void EffectiveResistanceDistance::firstJoin(std::vector<node> vecOfNodes,node v, node w,double edgeWeight){
    node u;
    for(count i=0;i<vecOfNodes.size();i++){
      u=vecOfNodes[i];
      M(u,v)=M(u,w)+edgeWeight;
      M(v,u)=M(u,v);
    }
  }
  void EffectiveResistanceDistance::edgeFire(std::vector<node> vecOfNodes,node v, node w,double edgeWeight){
    node x;
    node y;
    double fixFactor;
    arma::Mat<double> M2;

    fixFactor=4.*(edgeWeight+M(v,w));
    M2 = M;

    for(count i=0;i<vecOfNodes.size();i++){
      x=vecOfNodes[i];
      for(count j=i+1;j<vecOfNodes.size();j++){
        y=vecOfNodes[j];
        M2(x,y)=M(x,w)-M(x,v);
        M2(x,y)-=M(w,y)-M(v,y);
        M2(x,y)*=M2(x,y);
        M2(x,y)/=fixFactor;
        M2(x,y)=M(x,y)-M2(x,y);
        M2(y,x)=M2(x,y);
      }
    }
    M=M2;
  }

/**************************************************************************/
  void EffectiveResistanceDistance::nonBridgeDelete(std::vector<node> vecOfNodes,node v, node w,double edgeWeight){
    node x;
    node y;
    double fixFactor;
    arma::Mat<double> M2;

    fixFactor=4.*(edgeWeight- M(v,w));
    M2 = M;

    for(count i=0;i<vecOfNodes.size();i++){
      x=vecOfNodes[i];
      for(count j=i+1;j<vecOfNodes.size();j++){
        y=vecOfNodes[j];
        M2(x,y)=M(x,w)-M(x,v);
        M2(x,y)-=M(w,y)-M(v,y);
        M2(x,y)*=M2(x,y);
        M2(x,y)/=fixFactor;
        M2(x,y)=M(x,y)+M2(x,y);
        M2(y,x)=M2(x,y);
      }
    }
    M=M2;
  }
  void EffectiveResistanceDistance::corollary(std::vector<node> vecOfNodes,node u, node v, node w,double edgeWeightOne,double edgeWeightTwo){
    node x;
    node y;
    double A;
    double B;
    double C;
    double D;
    double E;
    double F;
    double N;
    double FD;
    double AmC;
    arma::Mat<double> M2;

    B=M(w,v)-M(v,u)+M(w,u);
    D=4.*(edgeWeightOne+ M(w,v));
    E=4.*(edgeWeightTwo-M(w,u));
    N=D*(D*E+4*B*B);
    M2 = M;

    for(count i=0;i<vecOfNodes.size();i++){
      x=vecOfNodes[i];
      A=M(x,v)-M(x,w);
      F=M(x,u)-M(x,w);
      for(count j=i+1;j<vecOfNodes.size();j++){
        y=vecOfNodes[j];
        C=M(v,y)-M(w,y);
        F-=M(u,y)-M(w,y);
        FD=F*D;
        AmC=A-C;
        M2(x,y)=FD-2*(AmC)*B;
        M2(x,y)*=M2(x,y);
        M2(x,y)/=N;
        M2(x,y)-=AmC*AmC/D;
        M2(y,x)=M(x,y)+M2(x,y);
      }
    }
    M=M2;
  }
  void EffectiveResistanceDistance::uncoarseTriangle(std::vector<node> vecOfNodes,std::tuple<node,node,count,double,double,double> triangle){
    node c,s,u,w;
    double edgeWeightcs,edgeWeightcw,edgeWeightsw;
    arma::Mat<double> L;
    c=std::get<0>(triangle);
    s=std::get<1>(triangle);
    w=std::get<2>(triangle);
    edgeWeightcs=std::get<3>(triangle);
    edgeWeightcw=std::get<4>(triangle);
    edgeWeightsw=std::get<5>(triangle);
    edgeWeightsw-=1./(1./edgeWeightcs+1./edgeWeightcw);
    L.set_size(3,3);
    L.zeros();
    //s--c--w
    L(0,1)=-edgeWeightcs;
    L(1,0)=L(0,1);
    L(0,2)=-edgeWeightsw;
    L(2,0)=L(0,2);
    L(1,2)=-edgeWeightcw;
    L(2,1)=L(1,2);
    L(0,0)=edgeWeightcs+edgeWeightsw;
    L(1,1)=edgeWeightcs+edgeWeightcw;
    L(2,2)=edgeWeightcw+edgeWeightsw;

    L=arma::pinv(L);

    M(c,w)=L(1,1)+L(2,2)-2*L(1,2);
    M(w,c)=M(c,w);
    M(c,s)=L(0,0)+L(1,1)-2*L(0,1);
    M(s,c)=M(c,s);

    for(count i=0;i<vecOfNodes.size();i++){
      u=vecOfNodes[i];
      if(!(u==s)&&!(u==w)){
        M(u,c)=std::min(M(u,s)+M(s,c),M(u,w)+M(w,c));
        M(c,u)=M(u,c);
      }
    }
  }
} /* namespace NetworKit*/
