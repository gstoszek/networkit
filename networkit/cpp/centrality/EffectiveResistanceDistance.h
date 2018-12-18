/*
* EffectiveResistanceDistance.h
        *
        *      Author: gstoszek
*/

#ifndef EFFECTIVERESISTANCEDISTANCE_H_
#define EFFECTIVERESISTANCEDISTANCE_H_

#define ARMA_DONT_PRINT_ERRORS
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
       *@param Laplacian, expect a Laplacian of smaller or equal size as EffectiveResistanceDistanceMatrice
       *@param vecOfNodes, expect a vector of nodes with equal size as Laplacian
       */
      void computeFromLaplacian(std::vector<node> vecOfNodes,arma::Mat<double> Laplacian);
      /**
       *executes first join operation
       *
       *@param vecOfNodes, expect a vector of nodes with equal size as Laplacian
       *@param v existing
       *@param w new node
       *@param edgeWeight, weight of edge reintroduced
       */
      void firstJoin(std::vector<node> vecOfNodes,node v, node w,double edgeWeight);
      /**
       *executes edge fire operation
       *
       *@param vecOfNodes, expect a vector of nodes with equal size as Laplacian
       *@param v and w, nodes of added edges
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
       *executes edge fire and non bridge delete operation in one step
       *
       *@param vecOfNodes, expect a vector of nodes with equal size as Laplacian
       *@param u,v and w, nodes of added edge [v,w] and deleted edge [u,w]
       *@param edgeWeightOne, weight of edge [v,w]
       *@param edgeWeightTwo, weight of edge [u,w]
       */
      void corollary(std::vector<node> vecOfNodes,node u, node v, node w,double edgeWeightOne,double edgeWeightTwo);
      /**
      *Executes uncoarsening of a inner innerTriangle;
      *@param c,s,w triangle nodes,
      *@param edgeWeightcs,edgeWeightcw,edgeWeightsw corresponding Weights.
      **/
      void uncoarseTriangle(std::vector<node> vecOfNodes,std::tuple<node,node,count,double,double,double> triangle);
      /**
       *EffectiveResistanceDistanceMatrice
       */
      arma::Mat<double> M;
    };


} /* namespace NetworKit */
#endif /* EFFECTIVERESISTANCEDISTANCE_H_*/
