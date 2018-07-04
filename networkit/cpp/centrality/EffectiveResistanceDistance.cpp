/*
*  EffectiveResistanceDistance.cpp
*
*      Author: gstoszek
*/

#include "EffectiveResistanceDistance.h"
#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include "../auxiliary/Log.h"
#include <chrono>
#include <stdlib.h>
#include <cmath>

namespace NetworKit {
/*
 * ***SAMPLING*** (Number *** II ***)
 */

  EffectiveResistanceDistance::EffectiveResistanceDistance(const Graph &G): Centrality(G, true) {

        D.resize(G.upperNodeIdBound());
        G.forNodes([&](node v){
            D[v].resize(G.upperNodeIdBound(), -1.);
            D[v][v]=0.;
        });
    }
    /*******************************************************************************************************************************************************/
    /*******************************************************************************************************************************************************/
    std::vector<std::vector<node>> EffectiveResistanceDistance::divide_and_conquer_distances(std::vector<std::vector<node>> List) {
        if (List.size() == 1) {
            for (count i = 1; i < List[0].size(); i++) {
                D[List[0][0]][List[0][i]] = 1.;
                D[List[0][i]][List[0][0]] = 1.;
            };
            for (count i = 1; i < List[0].size(); i++) {
                for (count j = i + 1; j < List[0].size(); j++) {
                    D[List[0][i]][List[0][j]] = 2.;
                    D[List[0][j]][List[0][i]] = 2.;
                }
            }
            return List;
        }
        else {
            std::vector<std::vector<node>> LeftList;
            std::vector<std::vector<node>> RightList;

            LeftList.resize(List.size()*0.5);
            RightList.resize(List.size()-LeftList.size());

            for(count i=0;i<LeftList.size();i++){
                LeftList[i]=List[i];
            }
            for(count i=0;i<RightList.size();i++){
                RightList[i]=List[LeftList.size()+i];
            }

            LeftList=divide_and_conquer_distances(LeftList);
            RightList=divide_and_conquer_distances(RightList);

            std::cout << "\n";
            for(count x=0;x<D.size();x++){
                for(count y=0;y<D.size();y++){
                    std::cout << D[x][y] <<"  ";
                }
                std::cout << "\n";
            }
            std::cout << "\n";

            return merge(LeftList,RightList);
        }
    }
    /*******************************************************************************************************************************************************/
    /*******************************************************************************************************************************************************/
    std::vector<std::vector<node>> EffectiveResistanceDistance::merge(std::vector<std::vector<node>> LeftList, std::vector<std::vector<node>> RightList) {
        std::vector<std::vector<node>> NewList;
        NewList.resize(LeftList.size()+RightList.size());
        bool first_join_found =false;

        /*Left Matrix*/
        for(count L_i=0;L_i<LeftList.size();L_i++) {
            for (count L_j = 0; L_j < LeftList[L_i].size(); L_j++) {
            /*Right Matrix*/
                for (count R_i = 0; R_i < RightList.size(); R_i++) {
                    for (count R_j = 0; R_j < RightList[R_i].size(); R_j++) {
                        /*check if connection between two Moduls exists*/
                        if (G.hasEdge(LeftList[L_i][L_j], RightList[R_i][R_j])) {
                            if (!first_join_found) {
                                first_join(LeftList[L_i][L_j], RightList[R_i][R_j],LeftList, RightList);
                                /*Link the two matrices together*/
                                for(count L_i=0;L_i<LeftList.size();L_i++) {
                                    NewList[L_i]=LeftList[L_i];
                                }
                                for(count R_i=0;R_i<RightList.size();R_i++){
                                    NewList[LeftList.size()+R_i]=RightList[R_i];
                                }
                                first_join_found = true;
                            } else {
                                edge_fire(LeftList[L_i][L_j], RightList[R_i][R_j], NewList);
                            }
                        }
                    }
                }
            }
        }
        /*Add Edge in case there is no edge beetween LeftList and RightList*/
        if(!first_join_found){
            first_join(LeftList[0][0],RightList[0][0],LeftList,RightList);
            std::pair <node,node> new_edge;
            new_edge.first = LeftList[0][0];
            new_edge.second= RightList[0][0];
            List_new_Edges.push_back(new_edge);
            for(count L_i=0;L_i<LeftList.size();L_i++) {
                NewList[L_i]=LeftList[L_i];
            }
            for(count R_i=0;R_i<RightList.size();R_i++){
                NewList[LeftList.size()+R_i]=RightList[R_i];
            }
        }
        return NewList;
    }
    /*******************************************************************************************************************************************************/
    /*******************************************************************************************************************************************************/
    void EffectiveResistanceDistance::first_join(node i, node j, std::vector<std::vector<node>> LeftList, std::vector<std::vector<node>> RightList) {
        for (count L_x = 0; L_x < LeftList.size(); L_x++) {
            for (count L_y = 0; L_y < LeftList[L_x].size(); L_y++) {
                for (count R_x = 0; R_x < RightList.size(); R_x++) {
                    for (count R_y = 0; R_y < RightList[R_x].size(); R_y++) {
                        D[LeftList[L_x][L_y]][RightList[R_x][R_y]] =
                                D[LeftList[L_x][L_y]][i] + 1 + D[j][RightList[R_x][R_y]];
                        D[RightList[R_x][R_y]][LeftList[L_x][L_y]] = D[LeftList[L_x][L_y]][RightList[R_x][R_y]];
                    }
                }
            }
        }
        std::cout << "Matrix D after first join of" << "\n";
        std::cout << "LeftList:" << "\n";
        for (count x = 0; x < LeftList.size(); x++) {
            for(count y=0; y<LeftList[x].size();y++){
                std::cout << LeftList[x][y] << "  ";
            }
            std::cout << "\n";
        }
        std::cout << "RightList:" << "\n";
        for (count x = 0; x < RightList.size(); x++) {
            for(count y=0; y<RightList[x].size();y++){
                std::cout << RightList[x][y] << "  ";
            }
            std::cout << "\n";
        }
        std::cout <<"Edge:[" << i <<","<< j << "]"<<"\n";
        for(count x=0;x<D.size();x++){
            for(count y=0;y<D.size();y++){
                std::cout << D[x][y] <<"  ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    /*******************************************************************************************************************************************************/
    /*******************************************************************************************************************************************************/
    void EffectiveResistanceDistance::edge_fire(node i, node j, std::vector<std::vector<node>> List) {
        double edgefire;
        std::vector<std::vector<double>> D2;
        int D_size;
        count D2_x;
        count D2_y;
        D_size=0;
        for(count x=0;x<List.size();x++) {
            D_size += List[x].size();
        }

        D2.resize(D_size);
        for(count x=0;x<D2.size();x++){
            D2[x].resize(x);
        }
        D2_x=0;
        D2_y=0;
        for(count x=0;x<List.size();x++) {
            for (count y = 0; y < List[x].size(); y++) {
                if(!((x==List.size()-1)&&(y==List[x].size()-1))) {
                    D2_x = D2_y + 1;
                    for (count x_2 = x; x_2 < List.size(); x_2++) {
                        for (count y_2 = 0; y_2 < List[x_2].size(); y_2++) {
                            if (!((x == x_2) && (y_2 <= y))) {
                                edgefire = (D[List[x][y]][j] - D[List[x][y]][i]) -
                                           (D[j][List[x_2][y_2]] - D[i][List[x_2][y_2]]);
                                edgefire = edgefire * edgefire / (4 + 4 * D[i][j]);
                                D2[D2_x][D2_y] = D[List[x][y]][List[x_2][y_2]] - edgefire;
                                D2_x++;
                                if((List[x][y]==1)&&(List[x_2][y_2]==2)){
                                    std::cout << "edgefire="<<edgefire<<"\n";
                                }
                            }
                        }
                    }
                }
                D2_y++;
            }
        }
        D2_x=0;
        D2_y=0;
        for(count x=0; x<List.size();x++) {
            for (count y = 0; y < List[x].size(); y++) {
                if(!((x==List.size()-1)&&(y==List[x].size()-1))) {
                    D2_x = D2_y + 1;
                    for (count x_2 = x; x_2 < List.size(); x_2++) {
                        for (count y_2 = 0; y_2 < List[x_2].size(); y_2++) {
                            if (!((x == x_2) && (y_2 <= y))) {
                                D[List[x][y]][List[x_2][y_2]] = D2[D2_x][D2_y];
                                D[List[x_2][y_2]][List[x][y]] = D2[D2_x][D2_y];
                                D2_x++;
                            }
                        }
                    }
                    D2_y++;
                }
            }
        }
    }
    /*******************************************************************************************************************************************************/
    /*******************************************************************************************************************************************************/
    void EffectiveResistanceDistance::clean_edges() {
        std::vector<node> Edge_to_clean;
        std::vector<std::vector<double>> D2;
        D2=D;
        double Non_bridge_delete;
        for(count i=0; i<List_new_Edges.size();i++) {
            for (count x = 0; x < D.size(); x++) {
                for (count y = x + 1; y < D.size(); y++) {
                    Non_bridge_delete = (D[x][List_new_Edges[i].second - 1] - D[x][List_new_Edges[i].first - 1]) - (D[List_new_Edges[i].second - 1][y] - D[List_new_Edges[i].first - 1][y]);
                    Non_bridge_delete = Non_bridge_delete * Non_bridge_delete / (4 - 4 * D[List_new_Edges[i].first - 1][List_new_Edges[i].second - 1]);
                    D2[x][y] = D[x][y] + Non_bridge_delete;
                    D2[y][x] = D2[x][y];
                }
            }
            D = D2;
        }
    }
    /*******************************************************************************************************************************************************/
    /*******************************************************************************************************************************************************/
    std::vector<std::vector<double>> EffectiveResistanceDistance::getEffectiveResistanceDistanceMatrix() {
        return D;
    }
    /*******************************************************************************************************************************************************/
    /*******************************************************************************************************************************************************/
   void EffectiveResistanceDistance::run() {

        std::vector<std::vector<node>> List;

        List.resize((count)(log(G.upperNodeIdBound())));
        for(count i=0; i<List.size();i++){
            List[i].resize(1);
        }

        std::vector<node> ResidualList;
        std::vector<std::vector<node>> Residual;
        Residual.resize(1);
        Residual[0].resize(1);
        ResidualList.resize(G.upperNodeIdBound());
        node max_node;
        double max_degree;
        int max_node_index;
        count k;

        G.forNodes([&](node v){
           ResidualList[v]=v;
        });
        /*Prepering List for Divide and Concquer*/
        for(count i=0; i<(count)(log(G.upperNodeIdBound()));i++) {
            max_node = 0;
            max_degree = 0.;
            for (count j = 0; j < ResidualList.size(); j++) {
                if (G.weightedDegree(ResidualList[j]) > max_degree) {
                    max_degree = G.weightedDegree(ResidualList[j]+1);
                    max_node = ResidualList[j];
                    max_node_index= (int) j;
                }
            }
            List[i][0]=max_node;
            ResidualList.erase (ResidualList.begin()+max_node_index);

            k=0;
            while(k<ResidualList.size()){
                if(G.hasEdge(List[i][0],ResidualList[k])) {
                    List[i].push_back (ResidualList[k]);
                    ResidualList.erase (ResidualList.begin()+k);
                }
                else{
                    k++;
                }
            }
        }

        for(count x=0;x<List.size();x++){
            for(count y=0;y<List[x].size();y++){
                std::cout << List[x][y] <<"  ";

            }
            std::cout << "\n";
        }

        /*Starting Divide and Conquer*/
        divide_and_conquer_distances(List);
        /*Connecting Residals*/
        for (count i=0;i<ResidualList.size();i++){
            Residual[0][0]=ResidualList[i];
            List=merge(List,Residual);
        }
        /*Cleaning temporary edges*/
        std::cout << "List_new_Edges.size()="<< List_new_Edges.size() << "\n";
        clean_edges();
    }
} /* namespace NetworKit*/