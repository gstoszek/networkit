/*
* ERDLEVEL.h
*
*      Author: gstoszek
*/

#ifndef ERDLEVEL_H_
#define ERDLEVEL_H_

#include "Centrality.h"
#include "ERD2.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include <armadillo>

namespace NetworKit {

    class ERDLevel
    {
    public:

        ERDLevel(count ID, std::vector<node> vList,std::vector<std::pair<node,node>> Matching, arma::Mat<double> Adj);

        void set_ID(count ID);
        void set_vList(std::vector<node> vList);
        void set_Matching(std::vector<std::pair<node,node>> Matching);
        void set_Adj(arma::Mat<double> Adj);

        void print();
        count get_ID();
        std::vector<count> get_vList();
        std::vector<std::pair<node,node>> get_Matching();
        arma::Mat<double> get_Adj();
        arma::vec get_CalvAdj(node v);

        ~ERDLevel() {};

    private:

        count LevelID;
        std::vector<node> LevelvList;
        std::vector<std::pair<count,count>> LevelMatching;
        arma::Mat<double> LevelAdj;
    };


} /* namespace NetworKit */
#endif /* ERDLEVEL_H_*/
