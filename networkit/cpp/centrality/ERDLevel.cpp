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

  ERDLevel::ERDLevel(count ID,std::vector<std::tuple<node,node,node,double,double,double>> vecOfTriangles){
    LevelID=ID;
    LevelVecOfTriangles=vecOfTriangles;
  }
  void ERDLevel::set(count ID,std::vector<std::tuple<node,node,node,double,double,double>> vecOfTriangles){
    LevelID=ID;
    LevelVecOfTriangles=vecOfTriangles;
  }
  void ERDLevel::setID(count ID){
    LevelID=ID;
  }
  void ERDLevel::setVecOfTriangles(std::vector<std::tuple<node,node,node,double,double,double>> vecOfTriangles){
    LevelVecOfTriangles=vecOfTriangles;
  }
  count ERDLevel::getID(){
    return LevelID;
  }
  std::vector<std::tuple<node,node,node,double,double,double>> ERDLevel::getVecOfTriangles(){
    return LevelVecOfTriangles;
  }
} /* namespace NetworKit*/
