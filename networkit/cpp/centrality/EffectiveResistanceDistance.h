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

namespace NetworKit {
    /*
     *
     */
    class EffectiveResistanceDistance: public NetworKit::Centrality {

    public:

        EffectiveResistanceDistance(const Graph& G);

        /*updated*/
        void run();
        std::vector<std::vector<node>> divide_and_conquer_distances(std::vector<std::vector<node>> List);
        std::vector<std::vector<node>> merge(std::vector<std::vector<node>> LeftList,std::vector<std::vector<node>> RightList);
        void first_join(node i, node j, std::vector<std::vector<node>> LeftList, std::vector<std::vector<node>> RightList);
        void edge_fire(node i,node j,std::vector<std::vector<node>> List);
        void clean_edges();
        std::vector<std::stack<node>> get_most_valuable_nodes(Graph& G);

        std::vector<std::vector<double>> getEffectiveResistanceDistanceMatrix();

    private:

        std::vector<std::vector<double>> D;
        std::vector<std::pair<node,node>> List_new_Edges;

    };


} /* namespace NetworKit */
#endif /* EFFECTIVERESISTANCEDISTANCE_H_*/