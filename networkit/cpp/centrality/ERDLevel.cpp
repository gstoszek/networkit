/*
*  ERDLevel.cpp
*
*      Author: gstoszek
*/

#include "ERDLevel.h"
#include "ERD2.h"
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

  ERDLevel::ERDLevel(count ID, std::vector<node> vList,std::vector<std::pair<node,node>> Matching,std::vector<std::vector<node>> TopMatch, arma::Mat<double> Adj){
    LevelID=ID;
    LevelvList=vList;
    LevelMatching=Matching;
    LevelTopMatch=TopMatch;
    LevelAdj=Adj;
  }

  void ERDLevel::print(){
    std::cout << "ID= " << LevelID << "\n";
    std::cout <<  "Nodes: ";
    for(count i=0;i<LevelvList.size()-1;i++){
      std::cout <<  LevelvList[i] <<" , ";
    }
    std::cout <<  LevelvList[LevelvList.size()-1] <<"\n";

    if(LevelMatching.size()>0){
      for(count i=0;i<LevelMatching.size()-1;i++){
        std::cout <<  LevelMatching[i].first <<" , ";
      }
      std::cout <<  LevelMatching[LevelMatching.size()-1].first <<"\n";

      for(count i=0;i<LevelMatching.size()-1;i++){
        std::cout <<  LevelMatching[i].second <<" , ";
      }
      std::cout <<  LevelMatching[LevelMatching.size()-1].second <<"\n";
    }
    LevelAdj.print("Adj:");
  }
  void ERDLevel::set(count ID, std::vector<node> vList,std::vector<std::pair<node,node>> Matching,std::vector<std::vector<node>> TopMatch, arma::Mat<double> Adj){
    LevelID=ID;
    LevelvList=vList;
    LevelMatching=Matching;
    LevelTopMatch=TopMatch;
    LevelAdj=Adj;
  }
  void ERDLevel::set_ID(count ID){
    LevelID=ID;
  }

  void ERDLevel::set_vList(std::vector<node> vList){
    LevelvList=vList;
  }

  void ERDLevel::set_Matching(std::vector<std::pair<node,node>> Matching){
    LevelMatching=Matching;
  }
  void ERDLevel::set_TopMatch(std::vector<std::vector<node>> TopMatch){
    LevelTopMatch=TopMatch;
  }

  void ERDLevel::set_Adj(arma::Mat<double> Adj){
    LevelAdj=Adj;
  }

  void ERDLevel::set_i_j_ofAdj(count i,count j,double value){
    LevelAdj(i,j)=value;
    LevelAdj(j,i)=value;
  }

  count ERDLevel::get_ID(){
    return LevelID;
  }
  std::vector<count> ERDLevel::get_vList(){
    return LevelvList;
  }
  std::vector<std::pair<node,node>> ERDLevel::get_Matching(){
    return LevelMatching;
  }
  std::vector<std::vector<node>> ERDLevel::get_TopMatch(){
    return LevelTopMatch;
  }
  std::vector<node> ERDLevel::get_vecofTopmatch(node v){
    return LevelTopMatch[v];
  }
  arma::Mat<double> ERDLevel::get_Adj(){
    return LevelAdj;
  }
  arma::vec ERDLevel::get_CalofAdj(node v){
    return LevelAdj.col(v);
  }
} /* namespace NetworKit*/
