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

        ERDLevel(count ID, std::vector<node> vecOfNodes,std::vector<std::tuple<node,node,count>> vecOfTransferredEdges);

        void set(count ID, std::vector<node> vecOfNodes,std::vector<std::tuple<node,node,count>> vecOfTransferredEdges);
        void setID(count ID);
        void setVecOfNodes(std::vector<node> vecOfNodes);
        void setVecOfTransferredEdges(std::vector<std::tuple<node,node,count>> vecOfTransferredEdges);

        count getID();
        std::vector<count> getVecOfNodes();
        std::vector<std::tuple<node,node,count>> getVecOfTransferredEdges();

        ~ERDLevel() {};

    private:

        count LevelID;
        std::vector<node> LevelVecOfNodes;
        std::vector<std::tuple<node,node,count>> LevelVecOfTransferredEdges;
    };


} /* namespace NetworKit */
#endif /* ERDLEVEL_H_*/
