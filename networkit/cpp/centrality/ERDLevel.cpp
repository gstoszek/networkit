/*
*  ERDLevel.cpp
*
*      Author: gstoszek
*/

#include "ERDLevel.h"
#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include "../auxiliary/Log.h"
#include <chrono>
#include <stdlib.h>
#include <cmath>
#include <armadillo>
#include <cstdlib>

namespace NetworKit {

  ERDLevel::ERDLevel(count ID, std::vector<node> vecOfNodes,std::vector<std::tuple<node,node,count>> vecOfTransferredEdges){
    LevelID=ID;
    LevelVecOfNodes=vecOfNodes;
    LevelVecOfTransferredEdges=vecOfTransferredEdges;
  }
  void ERDLevel::set(count ID, std::vector<node> vecOfNodes,std::vector<std::tuple<node,node,count>> vecOfTransferredEdges){
    LevelID=ID;
    LevelVecOfNodes=vecOfNodes;
    LevelVecOfTransferredEdges=vecOfTransferredEdges;
  }
  void ERDLevel::setID(count ID){
    LevelID=ID;
  }
  void ERDLevel::setVecOfNodes(std::vector<node> vecOfNodes){
    LevelVecOfNodes=vecOfNodes;
  }
  void ERDLevel::setVecOfTransferredEdges(std::vector<std::tuple<node,node,count>> vecOfTransferredEdges){
    LevelVecOfTransferredEdges=vecOfTransferredEdges;
  }
  count ERDLevel::getID(){
    return LevelID;
  }
  std::vector<count> ERDLevel::getVecOfNodes(){
    return LevelVecOfNodes;
  }
  std::vector<std::tuple<node,node,count>> ERDLevel::getVecOfTransferredEdges(){
    return LevelVecOfTransferredEdges;
  }
} /* namespace NetworKit*/
