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

namespace NetworKit {

    class ERDLevel
    {
    public:

        ERDLevel(node id,std::vector<node> vecOfNeighbors, std::vector<double> vecOfWeights);

        void set(node ID,std::vector<node> vecOfNeighbors, std::vector<double> vecOfWeights);
        void setNode(node ID);
        void setNeighbors(std::vector<node> vecOfNeighbours);
        void setWeights(std::vector<double> vecOfWeights);

        node id();
        std::vector<node> neighbors();
        std::vector<double> weights();

        ~ERDLevel() {};

    private:

        node ID;
        std::vector<node> nghbrs;
        std::vector<double> wghts;

    };


} /* namespace NetworKit */
#endif /* ERDLEVEL_H_*/
