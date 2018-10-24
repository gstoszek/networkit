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

        CurrentFlowGroupCloseness(const Graph& G, const double epsilon=0.1, const double delta=0.1, const count groupsize = 2,const count upperDegreeBound = 1);

        void run();

        void setCFGCC(double newCFGCC);
        std::vector<node> getNodesofGroup();
        double getCFGCC();
        void update_CFGCC(count n_peripheral_merges);

        void compute_initial_ERD(count upperDegreeBound);
        void compute_ERD();
        void first_join(node s, node v);
        void first_join_peripheral(node s, node v,count multiplier);
        void edge_fire(node v, node w);
        void non_bridge_delete(node v,node w);
        void coarse_L(std::vector<std::pair<count,count>> indices);
        void uncoarse_L();
        void clean_network();
        void uncoarse(node s,node v,count ID);

        std::vector<std::pair<node,node>> update_Matching(std::vector<std::pair<count,count>> indices);
        std::vector<std::pair<count,count>> peripheral_indices();
        std::vector<std::pair<count,count>> coarsing_indices(count cDegree,bool Random);
        std::vector<std::vector<node>> update_TopMatch(count minDegree);
        count update_minDegree(count minDegree);
        std::pair<node,node> random_edge();
        count merge_peripheral_nodes();


    private:
        count n;
        double epsilon;
        double delta;
        count groupsize;
        count upperDegreeBound;

        double CFGCC;
        std::vector<node> S;
        std::vector<node> vList;
        std::vector<std::vector<node>> TopMatch;
        std::vector<ERDLevel> LevelList;
        arma::Mat<double> ERD;
        arma::Mat<double> L;
        arma::Mat<double> Adj;
    };


} /* namespace NetworKit */
#endif /* CURRENTFLOWGROUPCLOSENESS_H_ */
