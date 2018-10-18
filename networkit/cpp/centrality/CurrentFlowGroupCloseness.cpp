/*
 * CurrentFlowGroupCloseness.cpp
 *
 *      Author: gstoszek
 */

#include "CurrentFlowGroupCloseness.h"
#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include "../auxiliary/Log.h"
#include <chrono>
#include <stdlib.h>
#include <cmath>
#include "EffectiveResistanceDistance.h"
#include "ERD2.h"
#include <armadillo>

namespace NetworKit {
/*
 * ***SAMPLING*** (Number *** II ***)
 */

   CurrentFlowGroupCloseness::CurrentFlowGroupCloseness(const Graph& G, const double epsilon, const double delta, const double universalConstant, const int groupsize) : Centrality(G, true), epsilon(epsilon), delta(delta), universalConstant(universalConstant), groupsize(groupsize){

        S.clear();
        S.resize(groupsize);

        D.resize(G.upperNodeIdBound()+1);
        G.forNodes([&](node v){
            D[v].resize(G.upperNodeIdBound()+1, -1.);
            D[v][v]=0.;
        });

        CFGCC = 0.;
    }
    /*******************************************************************************************************************************************************/
    /*******************************************************************************************************************************************************/
    std::vector<node> CurrentFlowGroupCloseness::getNodesofGroup(){
        return S;
    }
    /*******************************************************************************************************************************************************/
    /*******************************************************************************************************************************************************/
    double CurrentFlowGroupCloseness::getCFGCC() {
        return CFGCC;
    }
    /*******************************************************************************************************************************************************/
    /*******************************************************************************************************************************************************/
    void CurrentFlowGroupCloseness::setCFGCC(double newCFGCC) {
        CFGCC = newCFGCC;
    }
    /*******************************************************************************************************************************************************/
    /*******************************************************************************************************************************************************/

    void CurrentFlowGroupCloseness::run() {

        count n = G.upperNodeIdBound();

        //Vector p(n,0.0);
        //Vector b(n,0.0);

        std::vector<double> d;
        d.resize(n,(double) (n*n));
        Vector zeroVector(n, 0.0);
        std::vector<bool> T;
        double S_CFGCC=0.;
        double S_currentCFGCC;

        node s_next=0;


        std::cout << "Number of Nodes:" << n << "\n";

        Aux::Log::setLogLevel("DEBUG");

        /*******************************************************************************************************************************************************/
        /*Preparing Seed*/
        /*******************************************************************************************************************************************************/

        T.resize(n,true);
        /*
        count t_i;
        for (count i = 0; i < samplesize; i++) {
            t_i = rand() % n;
            while (std::find(Seed.begin(), Seed.end(), t_i) != Seed.end()) {
                t_i = rand() % n;
            }
            Seed[i] = t_i;
        }
        */
        /*
        for(count i=0;i<Seed.size();i++){
            Seed[i]=i;
        }
        */
        /*******************************************************************************************************************************************************/
        /*Solver Setup*/
        /*******************************************************************************************************************************************************/
        /*
        CSRMatrix L = CSRMatrix::laplacianMatrix(G);
        Lamg<CSRMatrix> Solver(1e-5);
        Solver.setupConnected(L);
        */
        /*******************************************************************************************************************************************************/
        /*calculating D-Matrix with Solver*/
        /******************************************************************************************************************************************************/
        /*for (count i= 0; i<n; i++) {
                b[i] = 1.;
                for (count j = 0; j < samplesize; j++) {
                    b[Seed[j]] = -1.;
                    p = zeroVector;
                    Solver.solve(b, p);
                    D[i][Seed[j]] = fabs(p[i] - p[Seed[j]]);
                    b[Seed[j]] = 0.;
            }
                b[i]=0.;
            }*/
        /**************************************************************************************************************/
        /*calculating D-Matrix with Divide and Conquer*/
        /**************************************************************************************************************/
        /*
        std::cout << "\n";
        std::cout << "Start: Computing ERD"<< "\n";

        EffectiveResistanceDistance ERD(G);

        auto start = std::chrono::high_resolution_clock::now();
        ERD.run();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end-start;
        std::cout << "Finished Calculation of ERD in " << diff.count() << "(s)" << "\n";

        D = ERD.getEffectiveResistanceDistanceMatrix();

        std::cout << "\n";
        */
        /**************************************************************************************************************/
        /*calculating D-Matrix with ERD2*/
        /**************************************************************************************************************/
        std::cout << "\n";
        std::cout << "Start: Computing ERD_2"<< "\n";
        ERD2 ERD(G);
        std::cout << "Constructor finished"<< "\n";
        auto start = std::chrono::high_resolution_clock::now();
        ERD.run();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end-start;
        std::cout << "Finished Calculation of ERD_2 in " << diff.count() << "(s)" << "\n";

        D = ERD.getERDMatrix();

        std::cout << "\n";

        /*Test-Prints*/
        /*
        std::cout << "\n";
        for(count x=0;x<D.size();x++){
            for(count y=0;y<D.size();y++){
                std::cout << D[x][y] <<"  ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
        */
        /**************************************************************************************************************/
        /*Greedy*/
        /**************************************************************************************************************/
        for(count i=0;i<S.size();i++){
            /*Maximal Gain Loop*/
            for (count j=0; j<T.size();j++) {
                //use vector of bools
                if (T[j]){
                    /*Sample Loop*/
                    S_currentCFGCC = 0.;
                    for (count l = 0; l < n; l++) {
                        if (D[l][j]< d[l])
                            S_currentCFGCC = S_currentCFGCC + D[l][j];
                        else {
                            S_currentCFGCC = S_currentCFGCC + d[l];
                        }
                    }
                    S_currentCFGCC = (double)(n) / S_currentCFGCC;
                    if (S_currentCFGCC > S_CFGCC) {
                        S_CFGCC = S_currentCFGCC;
                        s_next = j;
                    }
                }
            }
            for(count j= 0; j<T.size(); j++) {
                if (D[j][s_next] < d[j]) {
                    d[j] = D[j][s_next];
                }
            }
            S[i]=s_next+1;
            T[s_next]=false;
        }
        CFGCC=S_CFGCC;
    }


} /* namespace NetworKit*/
