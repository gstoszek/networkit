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

        ERDLevel(count ID, std::vector<node> vList,std::vector<std::pair<node,node>> Matching,std::vector<std::vector<node>> TopMatch, arma::Mat<double> Adj);

        void set(count ID, std::vector<node> vList,std::vector<std::pair<node,node>> Matching,std::vector<std::vector<node>> TopMatch, arma::Mat<double> Adj);
        void set_ID(count ID);
        void set_vList(std::vector<node> vList);
        void set_Matching(std::vector<std::pair<node,node>> Matching);
        void set_TopMatch(std::vector<std::vector<node>> TopMatch);
        void set_Adj(arma::Mat<double> Adj);
        void set_i_j_ofAdj(count i,count j,double value);

        void print();
        count get_ID();
        std::vector<count> get_vList();
        std::vector<std::pair<node,node>> get_Matching();
        std::vector<std::vector<node>> get_TopMatch();
        std::vector<node> get_vecofTopmatch(node v);
        arma::Mat<double> get_Adj();
        arma::vec get_CalofAdj(node v);

        ~ERDLevel() {};

    private:

        count LevelID;
        std::vector<node> LevelvList;
        std::vector<std::pair<node,node>> LevelMatching;
        std::vector<std::vector<node>> LevelTopMatch;
        arma::Mat<double> LevelAdj;
    };


} /* namespace NetworKit */
#endif /* ERDLEVEL_H_*/
