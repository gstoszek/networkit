/*
 * CurrentFlowGroupCloseness.h
 *
 *  Created on:
 *      Author: GSTOSZEK
 */

#include "CurrentFlowGroupCloseness.h"
#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"


/**
 * CSRMatrix.h for Laplacian
 */

namespace NetworKit {
/*
 * ***SAMPLING*** (Number *** II ***)
 */

   CurrentFlowGroupCloseness::CurrentFlowGroupCloseness(const Graph& G, const double epsilon, const double delta, const double universalConstant, const int groupsize) : Centrality(G, true), epsilon(epsilon), delta(delta), universalConstant(universalConstant), groupsize(groupsize){

    }



    void CurrentFlowGroupCloseness::run() {
        //GREEDY WITH SAMPLING
        count n = G.upperNodeIdBound();

        std::vector<node> stack_T(n);
        std::vector<node> stack_S(n);
        std::vector<double> p(n);
        std::vector<double> b(n);
        d.clear();
        d.resize(n, 0);

        double currentCFCC=0;
        double nextCFCC=0;

        node nextnode;
        count nodetoerase_index;
        count t_i;
        count samplesize= ceil(universalConstant * epsilon * epsilon * floor(log(n)));

        CSRMatrix L = CSRMatrix::laplacianMatrix(G);

        //Group Loop
        for (count i = 1; i <= groupsize; i++){
            //Maximal Gain Loop

            for (count j= 0; j<stack_T.size(); j++) {

                //main function
                if (std::find(stack_S.begin(), stack_S.end(), stack_T.at(j)) == stack_S.end()) {
                    //Sample
                    for (count l = 0; l < samplesize; l++) {
                        t_i = rand() % (n + 1);  //delete?
                        do {
                            t_i = rand() % (n + 1);
                        } while ((t_i == j) &&
                                 (std::find(stack_S.begin(), stack_S.end(), stack_T.at(t_i)) != stack_S.end()));

                        for(count a=0;a<n;a++){
                            b.at(a)=0.;
                            p.at(a)=0.;
                        }
                        b.at(j)=1.;
                        b.at(t_i)=-1.;

                        Lamg::setupConnected(L);
                        Lamg::parallelSolve(b,p);

                        if(p.at(j)-p.at(t_i)<d.at(t_i))
                            currentCFCC = currentCFCC + p.at(j)-p.at(t_i);
                        else{
                            currentCFCC = currentCFCC + d.at(t_i);
                        }
                    }
                    currentCFCC = samplesize/n*(stack_T.size())/currentCFCC;
                }
                    if (currentCFCC > nextCFCC) {
                        nextCFCC = currentCFCC;
                        currentCFCC = 0;
                        nextnode = stack_T.at(j);
                        nodetoerase_index = j;
                        }
                }
            stack_S.push_back(nextnode);
            //stack_T.erase(nodetoerase_index);

            //Update d


        }

    }



} /* namespace NetworKit*/