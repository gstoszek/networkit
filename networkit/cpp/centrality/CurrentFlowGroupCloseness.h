/*
 * CurrentFlowGroupCloseness.h
 *
 *      Author: gstoszek
 */

#ifndef CURRENTFLOWGROUPCLOSENESS_H_
#define CURRENTFLOWGROUPCLOSENESS_H_

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
    class CurrentFlowGroupCloseness: public NetworKit::Centrality {

    public:
        /**
         *
         * @param G An connected unweighted graph.
         * @param k Size of the group of nodes
         * @param CB If equal 0 runs simply algorithm without coursing, atherwise sets a Coarsening Bound
         */
        CurrentFlowGroupCloseness(const Graph& G,const count k = 2,const count CB = 2,const double epsilon=0.1);
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

        count k=2;
        count CB=2;
        count n;

        double CFGCC;
        double epsilon;
        std::vector<node> S;

        std::vector<node> vecOfNodes;
        std::vector<std::vector<node>> TopMatch;
        std::vector<ERDLevel> LevelList;

        EffectiveResistanceDistance ERD;

        arma::Mat<double> L;

        void greedy(count n_peripheral_merges);
        std::vector<std::vector<node>> updateTopMatch(count minDegree);
        count updateMinDegree(count minDegree);
        std::vector<std::pair<node,node>> updateMatching(std::vector<std::pair<count,count>> indices);
        std::vector<std::pair<count,count>> peripheralCoarsingIndices();
        std::vector<std::pair<count,count>> coarsingIndices(count cDegree,bool Random);
        void uncoarseEfffectiveResistanceDistanceMatrix(node s,node v,count ID);
        count mergePeripheralNodes();
        void coarseLaplacian(std::vector<std::pair<count,count>> indices);


    };

} /* namespace NetworKit */
#endif /* CURRENTFLOWGROUPCLOSENESS_H_ */
