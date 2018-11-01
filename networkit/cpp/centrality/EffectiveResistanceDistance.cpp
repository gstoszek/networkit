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

  void EffectiveResistanceDistance::computeFromPinvL(arma::Mat<double> PinvLaplacian,std::vector<node> vecOfNodes){
    node v;
    node w;
    /*calculate initial M Matrix*/
    for(count i=0;i<vecOfNodes.size();i++){
      v=vecOfNodes[i];
      for(count j=i+1;j<vecOfNodes.size();j++){
        w=vecOfNodes[j];
        M(v,w)=PinvLaplacian(i,i)+PinvLaplacian(j,j)-2.*PinvLaplacian(i,j);
        M(w,v)=M(v,w);
      }
    }
  }
  void EffectiveResistanceDistance::firstJoin(std::vector<node> vecOfNodes,node v, node w,double edgeWeight){
    node u;
    for(count i=0;i<vecOfNodes.size();i++){
      u=vecOfNodes[i];
      M(u,w)=M(u,v)+edgeWeight;
      M(w,u)=M(u,w);
    }
  }
  void EffectiveResistanceDistance::edgeFire(std::vector<node> vecOfNodes,node v, node w,double edgeWeight){
    node x;
    node y;
    double fixFactor;
    arma::Mat<double> M2;

    fixFactor=4.*(1.+M(v,w));
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

    fixFactor=4.*(1.- M(v,w));
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
      FD=F*D;
      for(count j=i+1;j<vecOfNodes.size();j++){
        y=vecOfNodes[j];
        C=M(u,y)-M(w,y);
        F-=M(u,y)-M(w,y);
        AmC=(A-C);
        M2(x,y)=FD-2*(AmC)*B;
        M2(x,y)*=M2(x,y);
        M2(x,y)/=N;
        M2(x,y)-=AmC*AmC/D;
        M2(y,x)=M(x,y)+M2(x,y);
      }
    }
    M=M2;
  }
} /* namespace NetworKit*/
