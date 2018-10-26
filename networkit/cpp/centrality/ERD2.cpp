/*
*  ERD2.cpp
*
*      Author: gstoszek
*/

#include "ERD2.h"
#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include "../auxiliary/Log.h"
#include <chrono>
#include <stdlib.h>
#include <cmath>
#include <armadillo>

namespace NetworKit {

  ERD2::ERD2(const Graph &G): Centrality(G, true) {
        n=G.upperNodeIdBound();
        /* EffectiveResistanceDistance-Matrix*/
        ERD.resize(n);
        /*Laplacian*/
        L.set_size(n,n);
        L.zeros();
        /* */
        vList.resize(n);
        vAdj_List.resize(n);

        G.forNodes([&](node v){
            ERD[v].resize(n);
            ERD[v][v]=0.;
            vAdj_List[v].resize(0);
            vList[v]=v;
            G.forNodes([&](node w){
              if(G.hasEdge(v,w)){
                /*instead of 1 one should use a function for weighted edges*/
                L(v,w)=-1.;
                L(v,v)+=1.;
                vAdj_List[v].push_back(w);
              }
            });
            nEdges+=L(v,v);
        });
      }
    /**************************************************************************/
    std::vector<std::vector<double>> ERD2::getERDMatrix() {
        return ERD;
    }
    /**************************************************************************/
    void ERD2::run(){
      /**Set LEVEL here**/
      upperLevelIdBound =1;
      count current_Level;
      current_Level=0;
      coarsed_List.resize(0);
      auto start = std::chrono::high_resolution_clock::now();
      while(current_Level<upperLevelIdBound){
        coarse_L(current_Level);
        current_Level+=1;
      }
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = end-start;
      std::cout<< "Coarsening finished in: " << diff.count() << "(s)" << "\n";
      /*invert*/
      start = std::chrono::high_resolution_clock::now();
      L=arma::pinv(L, 0.01);
      end = std::chrono::high_resolution_clock::now();
      diff = end-start;
      std::cout<< "Pinv() calculated in: " << diff.count() << "(s)" << "\n";
      /*calculate initial ERD Matrix*/
      for(count i=0;i<vList.size();i++){
        for(count j=i+1;j<vList.size();j++){
          ERD[vList[i]][vList[j]]=L(i,i)+L(j,j)-2.*L(i,j);
          ERD[vList[j]][vList[i]]=ERD[vList[i]][vList[j]];
        }
      }
      for(count i=0;i<coarsed_List.size();i++){
        for(count j=0;j<ERD.size();j++){
          ERD[coarsed_List[i].first][j]=ERD[vAdj_List[coarsed_List[i].first][0]][j]+1.;
          ERD[j][coarsed_List[i].first]=ERD[coarsed_List[i].first][j];
        }
      }
    }
    /***************************************************************************/
    void uncoarse_L(){

    }


    /***************************************************************************/
    void ERD2::coarse_L(count current_Level){
      std::vector<node> coarse_List;
      coarse_List.resize(0);
      for(count i=0;i<L.n_rows;i++){
        if(L(i,i)==current_Level+1){
          coarse_List.push_back(i);
        }
      }
      for(count i=0;i<coarse_List.size();i++){
        for(count j=1;j<vAdj_List[coarse_List[i]].size();j++){
          if(L(vAdj_List[coarse_List[i]][0],vAdj_List[coarse_List[i]][j])!=-1){
            L(vAdj_List[coarse_List[i]][0],vAdj_List[coarse_List[i]][j])=-1.;
            L(vAdj_List[coarse_List[i]][j],vAdj_List[coarse_List[i]][0])=-1.;
            tM[vList[i]][vAdj_List[coarse_List[i]][j]];
            tM[vAdj_List[coarse_List[i]][j]][vList[i]];
            /*Diagonal*/
            L(vAdj_List[coarse_List[i]][0],vAdj_List[coarse_List[i]][0])+=1.;
          }
        }
        L(vAdj_List[coarse_List[i]][0],vAdj_List[coarse_List[i]][0])-=1.;
        coarsed_List.push_back(std::make_pair(vList[coarse_List[i]],current_Level));
      }
      arma::uvec indices(vList.size()-coarse_List.size());
      count k;
      count j;
      k=0;
      j=0;
      for(count i=0;i<vList.size();i++){
        if(coarse_List[k]==i){
          k++;
        }
        else{
          indices(j)=i;
          j++;
        }
      }
      for(count i=0;i<coarse_List.size();i++){
        vList.erase(vList.begin()+coarse_List[i]-i);
      }
      L=L.submat(indices, indices);
    }
    /**************************************************************************/
    /*
    void ERD2::sort_Level_List(std::vector<node> vLevel_List,count LevelId){
      std::vector<std::vector<node>> bucket_List;
      bucket_List.resize(upperLevelIdBound)
      bucket_List.resize(0);
      count k,
      k=0;
      while(k!=vLevel_List.size()){
        if(nAdj_List[vLevel_List[k]].second<upperNodeIdBound){
          bucket_List[nAdj_List[vLevel_List[k]].second-1].push_back vLevel_List[k];
        }
        else{
          bucket_List[upperNodeIdBound].push_back vLevel_List[k];
        }
        k++;
      }
      for(count i=0;i<bucket_List.size();i++){
        for(count j=0;j<bucket_List[i].size();j++){
          vLevel_List=bucket_List[i][j];
        }
      }
    }
    */
    /**************************************************************************/
    /*
    void ERD2::coarse(count samplesize, count LevelID){
      if(LevelID==3){
      }
      else{
        std::pair<node,node> coarsed_edge;
        std::vector<std::pair<bool,double>> Absorber_Edge_List;
        std::vector<std::pair<bool,double>> Submitter_Edge_List;
        Absorber_Edge_List.resize(Adj.size());
        Submitter_Edge_List.resize(Adj.size());
        coarsed_edge=select_edge();
        //coarsed_edge=random_edge();
        Absorber_Edge_List=Adj[coarsed_edge.first];
        Submitter_Edge_List=Adj[coarsed_edge.second];
        coarsening(coarsed_edge);
        coarse(1);
        uncoarsening(Absorber_Edge_List,Submitter_Edge_List,coarsed_edge);
      }
    }
    */
    /**************************************************************************/
    /*
    std::pair<node,node> ERD2::select_edge(){
      return make_pair(Matching_List[vLevel_List[0]].second,vLevel_List[0])
    }
    */
    /**************************************************************************/
    /*
    std::pair<node,node> ERD2::random_edge(){
      count Random_Edge;
      count i;
      count a;
      count b;

      a=1;
      b=nAdj_List.size()-1;

      i=a+(b-a)*0.5;
      Random_Edge= (std::rand() % nEdges)+1;
      while(!((Random_Edge>nAdj_List[i-1].second)&&(Random_Edge<=nAdj_List[i].second))){
        if(Random_Edge<=nAdj_List[i-1].second){
          b=i;
        }
        else{
          a=i;
        }
        i=a+(b-a)*0.5;
      }
      count k=0;
      count j=0;
      while(k<Random_Edge-nAdj_List[i-1].second){
        if(Adj[nAdj_List[i].first][j].first){
            k++;
        }
        j++;
      }
      j--;
      if(j>=G.upperNodeIdBound()){
        std::cout<<"<"<<nAdj_List[i].first<<","<<j<<"> mit "<<Random_Edge<<"\n";
        std::cout<<"No. of Edges="<<nEdges<<"\n";
        std::cout<<"i-1="<<i+1<<"-->"<<"v="<<nAdj_List[i-1].first<<" No.="<< nAdj_List[i-1].second<<"\n";
        std::cout<<"i="<<i<<"-->"<<"v="<<nAdj_List[i].first<<" No.="<< nAdj_List[i].second<<"\n";
        std::cout<<"i+1="<<i+1<<"-->"<<"v="<<nAdj_List[i+1].first<<" No.="<< nAdj_List[i+1].second<<"\n";
        count p;
        for(count k=0;k<Adj.size();k++){
          if(Adj[nAdj_List[i].first][k].first){
            p++;
          }
        }
        std::cout<<"p="<<p<<"\n";
      }
      return std::make_pair(nAdj_List[i].first,j);
    }
    */
    /**************************************************************************/
    /*
    void ERD2::coarsening_2(std::pair<node,node> coarsed_edge){
      std::cout << "coarsening! <a,s>=<" << coarsed_edge.first << "," << coarsed_edge.second <<">"<<"\n";
      vLevel_List.erase (vLevel_List.begin());
      Adj[coarsed_edge.first][coarsed_edge.second].first=false;
      Adj[coarsed_edge.second][coarsed_edge.first].first=false;
    }
    */
    /**************************************************************************/
    /**************************************************************************/
    /*void ERD2::coarsening(std::pair<node,node> coarsed_edge){
      count x;
      count y;
      std::vector<int> transfered_Edges;
      transfered_Edges.resize(Adj.size()+1);
      transfered_Edges[0]=0;
      for(count i=1;i<nAdj_List.size();i++){
        x=nAdj_List[i].first;
        if(Adj[coarsed_edge.second][x].first==true){
          if(Adj[coarsed_edge.first][x].first==false){
            Adj[coarsed_edge.first][x]=std::make_pair(true,Adj[coarsed_edge.second][x].second);
            Adj[x][coarsed_edge.first]=std::make_pair(true,Adj[coarsed_edge.second][x].second);
            transfered_Edges[x+1]++;
            transfered_Edges[coarsed_edge.first+1]++;
          }
          Adj[coarsed_edge.second][x]=std::make_pair(false,0.);
          Adj[x][coarsed_edge.second]=std::make_pair(false,0.);
          transfered_Edges[x+1]--;
          transfered_Edges[coarsed_edge.second+1]--;
        }
      }
      Adj[coarsed_edge.first][coarsed_edge.first]=std::make_pair(false,0);
      transfered_Edges[coarsed_edge.first+1]-=2;
      y=1;
      nAdj_List[1].second+=transfered_Edges[nAdj_List[1].first+1];
      for(count i=2;i<nAdj_List.size();i++){
          transfered_Edges[nAdj_List[i].first+1]+=transfered_Edges[nAdj_List[i-1].first+1];
          nAdj_List[i].second+=transfered_Edges[nAdj_List[i].first+1];
          if(nAdj_List[i].first==coarsed_edge.second){
            y=i;
          }
      }
      nAdj_List.erase(nAdj_List.begin()+y);
      nEdges=nAdj_List[nAdj_List.size()-1].second;
    }
    */
    /**************************************************************************/
    /*
    void ERD2::uncoarsening(std::vector<std::pair<bool,double>> Absorber_Edge_List,
      std::vector<std::pair<bool,double>> Submitter_Edge_List,
      std::pair<node,node> coarsed_edge){
        count x;

        first_join(coarsed_edge);
        nAdj_List.push_back(std::make_pair(coarsed_edge.second,1));
        for(count i=1; i<nAdj_List.size();i++){
          x=nAdj_List[i].first;
          if((Submitter_Edge_List[x].first) && (x!=coarsed_edge.first)){
            if(Absorber_Edge_List[x].first){
              edge_fire(std::make_pair(x,coarsed_edge.second));
            }
            else{
              ERD[coarsed_edge.second][x]=ERD[coarsed_edge.first]
            }
          }
          if((Adj[coarsed_edge.first][x].first)&&(Absorber_Edge_List[x].first!=true)){
            non_bridge_delete(std::make_pair(x,coarsed_edge.first));
          }
        }
      }
      */
    /**************************************************************************/
    /*
    void ERD2::first_join(std::pair<node,node> coarsed_edge){
      for(count i=1;i<nAdj_List.size();i++){
        ERD[nAdj_List[i].first][coarsed_edge.second]=ERD[nAdj_List[i].first][coarsed_edge.first]+1.;
        ERD[coarsed_edge.second][nAdj_List[i].first]=ERD[nAdj_List[i].first][coarsed_edge.second];
      }
      Adj[coarsed_edge.first][coarsed_edge.second]=std::make_pair(true,1.);
      Adj[coarsed_edge.second][coarsed_edge.first]=std::make_pair(true,1.);
    }
    */
    /**************************************************************************/
    /*
    void ERD2::edge_fire(std::pair<node,node> coarsed_edge){
      std::vector<std::vector<double>> ERD2;
      ERD2.resize(ERD.size());
      G.forNodes([&](node v){
        ERD2[v].resize(ERD.size());
      });
      count x;
      count y;

      for(count i=1;i<nAdj_List.size();i++){
        x=nAdj_List[i].first;
        for(count j=i+1;j<nAdj_List.size();j++){
          y=nAdj_List[j].first;
          ERD2[x][y]=ERD[x][coarsed_edge.second]-ERD[x][coarsed_edge.first];
          ERD2[x][y]-=ERD[coarsed_edge.second][y]-ERD[coarsed_edge.first][y];
          ERD2[x][y]*=ERD2[x][y];
          ERD2[x][y]/=4.*(1.+ERD[coarsed_edge.first][coarsed_edge.second]);
          ERD2[x][y]=ERD[x][y]-ERD2[x][y];
          ERD2[y][x]=ERD2[x][y];
        }
      }
      ERD=ERD2;
      Adj[coarsed_edge.first][coarsed_edge.second]=std::make_pair(true,1.);
      Adj[coarsed_edge.second][coarsed_edge.first]=std::make_pair(true,1.);
    }
    */
    /**************************************************************************/
    /*
    void ERD2::non_bridge_delete(std::pair<node,node> coarsed_edge){
      std::vector<std::vector<double>> ERD2;
      ERD2.resize(ERD.size());
      G.forNodes([&](node v){
        ERD2[v].resize(ERD.size());
      });
      count x;
      count y;

      for(count i=1;i<nAdj_List.size();i++){
        x=nAdj_List[i].first;
        for(count j=i+1;j<nAdj_List.size();j++){
          y=nAdj_List[j].first;
          ERD2[x][y]=ERD[x][coarsed_edge.second]-ERD[x][coarsed_edge.first];
          ERD2[x][y]-=ERD[coarsed_edge.second][y]-ERD[coarsed_edge.first][y];
          ERD2[x][y]*=ERD2[x][y];
          ERD2[x][y]/=4.*(1.-ERD[coarsed_edge.first][coarsed_edge.second]);
          ERD2[x][y]=ERD[x][y]+ERD2[x][y];
          ERD2[y][x]=ERD2[x][y];
        }
      }
      ERD=ERD2;
      Adj[coarsed_edge.first][coarsed_edge.second]=std::make_pair(false,0.);
      Adj[coarsed_edge.second][coarsed_edge.first]=std::make_pair(false,0.);
    }
    */
    /**************************************************************************/
} /* namespace NetworKit*/
