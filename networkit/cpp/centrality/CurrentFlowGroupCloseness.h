/*
 * CurrentFlowGroupCloseness.h
 *
 *      Author: gstoszek
 */

#ifndef CURRENTFLOWGROUPCLOSENESS_H_
#define CURRENTFLOWGROUPCLOSENESS_H_

#define ARMA_DONT_PRINT_ERRORS
#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include "../components/ConnectedComponents.h"
#include <armadillo>
#include <algorithm>

namespace NetworKit {
    /*
     *
     */
    class CurrentFlowGroupCloseness: public NetworKit::Algorithm {

    public:
        /**
         *
         * @param G An connected unweighted graph.
         * @param k Size of the group of nodes
         * @param CB If equal 0 runs simply algorithm without coursing, atherwise sets a Coarsening Bound
         */
        CurrentFlowGroupCloseness(Graph& G,const count k = 2,const count CB = 2,const double epsilon=0.1, const bool doInvert=true);
        /**
         * Computes group of size k with maximum closeness and coresponding value on the graph passed in the constructor.
         */
        void run();
        /**
         * Returns group of size k with maximum closeness.
         */
        std::vector<node> getNodesofGroup();
        /**
         * Returns maximum current flow group closeness
         */
        double getCFGCC();
        /**
         *Computes value for given group
         */



    private:

        Graph& G;
        count k=2;
        count CB=2;
        bool doInvert=true;

        double CFGCC;
        double epsilon;
        count n;
        count limit;

        std::vector<bool> E;
        std::vector<node> S;
        std::vector<node> coarsedNodes;
        std::vector<std::vector<node>> coarsedNeighbors;
        std::vector<std::vector<double>> coarsedWeights;
        std::vector<std::vector<double>> coarsedDistances;
        void greedy();
        void greedyLAMG();
        count minDegree();
        std::vector<node> coarsingNodes(count cDegree,bool Random);
        void graphCoarsening(std::vector<node> nodes);
        void adjacentCliqueWeights(std::vector<node> neighbors,std::vector<double> weights);
        arma::Mat<double> pinv();
        double exactDistance(count c,std::vector<double> dst, arma::Mat<double> Pinv);
        arma::Mat<double> distanceMatrix(std::vector<node> neighbors,arma::Mat<double> Pinv);
        double residual(double variation, double factor);
        double distanceFinder(node y, count x_i);
        double starWeight(count c);
    };

} /* namespace NetworKit */
#endif /* CURRENTFLOWGROUPCLOSENESS_H_ */
