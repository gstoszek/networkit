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
         */
        CurrentFlowGroupCloseness(Graph& G,const count k = 2,const double epsilon=0.1);
        /**
         * Computes group of size k with maximum closeness and coresponding value on the graph passed in the constructor.
         */
        void run();
        /**
         * Returns group of size k with maximum closeness.
         */
        std::vector<node> getNodesofGroup();
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

        double CFGCC;
        std::vector<node> S;

        void runExact();
        void runExactLS();
        /*Move to other class*/
        CSRMatrix shed(CSRMatrix matrix, count s);
        CSRMatrix gaussianMatrx(count numberOfRows,count numberOfColumns);
    };

} /* namespace NetworKit */
#endif /* CURRENTFLOWGROUPCLOSENESS_H_ */
