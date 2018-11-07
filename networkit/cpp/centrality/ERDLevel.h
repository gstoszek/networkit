/*
* ERDLEVEL.h
*
*      Author: gstoszek
*/

#ifndef ERDLEVEL_H_
#define ERDLEVEL_H_

#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include <armadillo>

namespace NetworKit {

    class ERDLevel
    {
    public:

        ERDLevel(count ID,std::vector<std::tuple<node,node,node,double,double,double>> vecOfTriangles);

        void set(count ID,std::vector<std::tuple<node,node,node,double,double,double>> vecOfTriangles);

        void setID(count ID);
        void setVecOfTriangles(std::vector<std::tuple<node,node,node,double,double,double>> vecOfTriangles);

        count getID();
        std::vector<std::tuple<node,node,node,double,double,double>> getVecOfTriangles();

        ~ERDLevel() {};

    private:

        count LevelID;
        std::vector<std::tuple<node,node,node,double,double,double>> LevelVecOfTriangles;
    };


} /* namespace NetworKit */
#endif /* ERDLEVEL_H_*/
