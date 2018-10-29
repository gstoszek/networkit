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

        ERDLevel(count ID, std::vector<node> vList,std::vector<std::pair<node,node>> Matching,std::vector<std::vector<node>> TopMatch);

        void set(count ID, std::vector<node> vList,std::vector<std::pair<node,node>> Matching,std::vector<std::vector<node>> TopMatch);
        void set_ID(count ID);
        void set_vList(std::vector<node> vList);
        void set_Matching(std::vector<std::pair<node,node>> Matching);
        void set_TopMatch(std::vector<std::vector<node>> TopMatch);

        void print();
        count get_ID();
        std::vector<count> get_vList();
        std::vector<std::pair<node,node>> get_Matching();
        std::vector<std::vector<node>> get_TopMatch();
        std::vector<node> get_vecofTopmatch(node v);

        ~ERDLevel() {};

    private:

        count LevelID;
        std::vector<node> LevelvList;
        std::vector<std::pair<node,node>> LevelMatching;
        std::vector<std::vector<node>> LevelTopMatch;
    };


} /* namespace NetworKit */
#endif /* ERDLEVEL_H_*/
