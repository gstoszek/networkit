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
#include "EffectiveResistanceDistance.h"
#include "ERDLevel.h"
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
        std::vector<node> vecOfPeripheralNodes;
        std::vector<node> S;
        std::vector<ERDLevel> coarsedNodes;
        void greedy();
        void greedyLAMG();
        count updateMinDegree();
        std::vector<node> coarsingIndices(count cDegree,bool Random);
        void uncoarseEfffectiveResistanceDistanceMatrix(count ID);
        void mergePeripheralNodes();
        void coarseGraph(std::vector<node> vecOfChosenNodes,count degree);
        void computeStarCliqueWeights(node c);
        arma::Mat<double> computePinvOfLaplacian();

    };

} /* namespace NetworKit */
#endif /* CURRENTFLOWGROUPCLOSENESS_H_ */
