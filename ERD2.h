/*
* ERD2.h
*
*      Author: gstoszek
*/

#ifndef ERD2_H_
#define ERD2_H_

#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include <armadillo>

namespace NetworKit {

    class ERD2: public NetworKit::Centrality {

    public:

        ERD2(const Graph& G);

        void run();
        std::vector<std::vector<double>> getERDMatrix();
        void coarsening(std::pair<node,node> coarsed_edge);
        void uncoarsening(std::vector<std::pair<bool,double>> Absorber_Edge_List,std::vector<std::pair<bool,double>> Submitter_Edge_List,std::pair<node,node> coarsed_edge);
        void first_join(std::pair<node,node> coarsed_edge);
        void edge_fire(std::pair<node,node> coarsed_edge);
        void non_bridge_delete(std::pair<node,node> coarsed_edge);
        void preprocessing();
        void coarse(count samplesize);
        std::pair<node,node> random_edge();
        void coarse_L(count current_Level);
        void uncoarse_L();

    private:
        std::vector<std::vector<double>> ERD;
        arma::Mat<double> L;
        std::vector<node> vList;
        std::vector<std::vector<node>> vAdj_List;
        std::vector<std::vector<bool>> tM;
        std::vector<std::pair<node,count>> coarsed_List;
        count nEdges;
        count upperLevelIdBound;
        count n;
    };


} /* namespace NetworKit */
#endif /* ERD2_H_*/
