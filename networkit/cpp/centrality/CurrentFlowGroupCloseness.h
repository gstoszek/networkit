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
#include "ERD2.h"
#include "ERDLevel.h"
#include <armadillo>

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
         * @
       */
        CurrentFlowGroupCloseness(const Graph& G,const count k = 2,const count CB = 2);
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
        std::vector<node> S;

        std::vector<node> vList;
        std::vector<std::vector<node>> TopMatch;
        std::vector<ERDLevel> LevelList;
        arma::Mat<double> ERD;
        arma::Mat<double> L;
        arma::Mat<double> Adj;

        void cleanNetwork();
        void greedy(count n_peripheral_merges);
        void computeInitialERD(count upperDegreeBound);
        void computeERD();
        std::vector<std::vector<node>> updateTopMatch(count minDegree);
        count updateMinDegree(count minDegree);
        std::vector<std::pair<node,node>> updateMatching(std::vector<std::pair<count,count>> indices);
        std::vector<std::pair<count,count>> peripheralCoarsingIndices();
        std::vector<std::pair<count,count>> coarsingIndices(count cDegree,bool Random);
        void uncoarse(node s,node v,count ID);
        count mergePeripheralNodes();
        void coarseLaplacian(std::vector<std::pair<count,count>> indices);
        void firstJoinPeripheral(node s, node v,count multiplier);
        void firstJoin(node s, node v);
        void edgeFire(node v, node w);
        void nonBridgeDelete(node v,node w);

    };

} /* namespace NetworKit */
#endif /* CURRENTFLOWGROUPCLOSENESS_H_ */
