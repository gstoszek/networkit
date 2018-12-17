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
        std::vector<node> vecOfPeripheralNodes;
        std::vector<node> S;
        std::vector<node> coarsedNodes;
        std::vector<std::vector<node>> coarsedNeighbors;
        std::vector<std::vector<double>> coarsedWeights;
        std::vector<std::vector<double>> coarsedDistances;
        void greedy();
        void greedyLAMG();
        count updateMinDegree();
        std::vector<node> coarsingIndices(count cDegree,bool Random);
        void mergePeripheralNodes();
        void coarseGraph(std::vector<node> vecOfChosenNodes,count degree);
        void computeStarCliqueWeights(node c);
        arma::Mat<double> computePinvOfLaplacian();
        double computeApproxDistances(std::vector<bool> W, std::vector<node> reverse, std::vector<double> dst, std::vector<double> *dstApprox);
        double computeExactDistance(std::vector<bool> E,std::vector<node> reverse, std::vector<double> dst, std::vector<double> *dstApprox, arma::Mat<double> Pinv);
        arma::Mat<double> computeDistanceMatrix(std::vector<node> neighbors,std::vector<count> neighbors_i,std::vector<bool> E,arma::Mat<double> Pinv);
        double residual(double dstxj, double dstxi, double dstjy,double dstiy,double dstij,double weight,double factor);
        double distanceFinder(node y, count x_i);
    };

} /* namespace NetworKit */
#endif /* CURRENTFLOWGROUPCLOSENESS_H_ */
