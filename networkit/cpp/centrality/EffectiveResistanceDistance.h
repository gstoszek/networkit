/*
* EffectiveResistanceDistance.h
        *
        *      Author: gstoszek
*/

#ifndef EFFECTIVERESISTANCEDISTANCE_H_
#define EFFECTIVERESISTANCEDISTANCE_H_

#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include <armadillo>

namespace NetworKit {

    class
    EffectiveResistanceDistance{

    public:
      /**
       *Constructor
       */
      EffectiveResistanceDistance(const count n = 0);
      /**
       *computes EffectiveResistanceDistanceMatrice
       *
       *@param Moore Penrose Pseudo Inverse Laplacian, expect a Laplacian of smaller or equal size as EffectiveResistanceDistanceMatrice
       *@param vecOfNodes, expect a vector of nodes with equal size as Laplacian
       */
      void computeFromPinvL(arma::Mat<double> PinvLaplacian,std::vector<node> vecOfNodes);
      /**
       *executes first join operation
       *
       *@param vecOfNodes, expect a vector of nodes with equal size as Laplacian
       *@param s and v, nodes of added edges
       *@param edgeWeight, weight of edge reintroduced
       */
      void firstJoin(std::vector<node> vecOfNodes,node v, node w,double edgeWeight);
      /**
       *executes edge fire operation
       *
       *@param vecOfNodes, expect a vector of nodes with equal size as Laplacian
       *@param s and v, nodes of added edges
       *@param edgeWeight, weight of edge reintroduced
       */
      void edgeFire(std::vector<node> vecOfNodes,node v, node w,double edgeWeight);
      /**
       *executes non bridge delete operation
       *
       *@param vecOfNodes, expect a vector of nodes with equal size as Laplacian
       *@param s and v, nodes of added edges
       *@param edgeWeight, weight of edge reintroduced
       */
      void nonBridgeDelete(std::vector<node> vecOfNodes,node v, node w,double edgeWeight);
      /**
       *EffectiveResistanceDistanceMatrice
       */
      arma::Mat<double> M;
    };


} /* namespace NetworKit */
#endif /* EFFECTIVERESISTANCEDISTANCE_H_*/
