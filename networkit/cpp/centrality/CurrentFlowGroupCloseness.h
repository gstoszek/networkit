/*
 * CurrentFlowGroupCloseness.h
 *
 *      Author: gstoszek
 */

#ifndef CURRENTFLOWGROUPCLOSENESS_H_
#define CURRENTFLOWGROUPCLOSENESS_H_

#define ARMA_DONT_PRINT_ERRORS
#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include "../numerics/ConjugateGradient.h"
#include "../components/ConnectedComponents.h"
#include <armadillo>
#include <algorithm>

namespace NetworKit {
    /*
     *
     */
    class CurrentFlowGroupCloseness: public NetworKit::Algorithm {

    public:
        /**
         *
         * @param G An connected unweighted graph.
         * @param k Size of the group of nodes
         * @param t type of computation (see below)
         */
        CurrentFlowGroupCloseness(Graph& G,const count k = 2,const double epsilon=0.1, const count t=1);
        /**
         * Computes group of size k with maximum closeness and coresponding value on the graph passed in the constructor.
         */
        void run();
        /**
         * Returns group of size k with maximum closeness.
         */
        std::vector<node> getNodesOfGroup();
        /**
         * Returns maximum current flow group closeness
         */
        double getCFGCC();
        /**
         *Computes value for given group
         */



    private:

        Graph& G;
        count k=2;
        double epsilon;
        count t=1;

        double CFGCC;
        std::vector<node> S;

        /*
        * t = 1, runs Determinisitc Greedy Algorithm for Current flow groupp closeness 1
        */
        void runDeterministicGreedy();
        /*
        * t = 2, runs Determinisitc Greedy Algorithm for Current flow groupp closeness 1,
        * which uses less space. !!!SLOW!!!
        */
        void runDeterministicGreedyLS();
        /*
        * t = 3, computes TopK with pseudo inverse for definition 1
        */
        void runTopK();
        /*
        * t = 4, runs approximation Greedy Algorithm for Current flow groupp closeness 1,
        * as proposed by Li et al. !!! Choose epsilon correctly!!!
        */
        void runApproxJLT();
        /*
        * t = 5, computes optimal solution
        */
        void runOptimum1();
        /*
        * uses t = 6 and t = 7,
        */
        void runRandom();
        /*
        * t = 6, computes closeness 1 of a given group s
        */
        void closeness1(std::vector<bool> s);
        /*
        * t = 7, computes closeness 2 of a given group s
        */
        void closeness2(std::vector<bool> s);
        /*
        * t = 8, runs Determinisitc Greedy Algorithm for Current flow groupp closeness 2
        */
        void runGCFCC2();

        void runApproxTopOne();
        /*Move to other class*/

        CSRMatrix shed(CSRMatrix matrix, count s);
        CSRMatrix gaussianMatrx(count numberOfRows,count numberOfColumns);
    };

} /* namespace NetworKit */
#endif /* CURRENTFLOWGROUPCLOSENESS_H_ */
