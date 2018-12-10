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
#include <cstdlib>

namespace NetworKit {

  ERDLevel::ERDLevel(node id,std::vector<node> vecOfNeighbors, std::vector<double> vecOfWeights){
    ID=id;
    nghbrs=vecOfNeighbors;
    wghts=vecOfWeights;
  }
  void ERDLevel::set(node id,std::vector<node> vecOfNeighbors, std::vector<double> vecOfWeights){
    ID=id;
    nghbrs=vecOfNeighbors;
    wghts=vecOfWeights;
  }
  void ERDLevel::setNode(node id){
    ID=id;
  }
  void ERDLevel::setNeighbors(std::vector<node> vecOfNeighbors){
    nghbrs=vecOfNeighbors;
  }
  void ERDLevel::setWeights(std::vector<double> vecOfWeights){
    wghts=vecOfWeights;
  }
  node ERDLevel::id(){
    return ID;
  }
  std::vector<node> ERDLevel::neighbors(){
    return nghbrs;
  }
  std::vector<double> ERDLevel::weights(){
    return wghts;
  }
} /* namespace NetworKit*/
