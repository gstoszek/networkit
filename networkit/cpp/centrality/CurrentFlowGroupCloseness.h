/*
 * CurrentFlowGroupCloseness.h
 *
 *  Created on:
 *      Author: GSTOSZEK
 */

#ifndef CURRENTFLOWGROUPCLOSENESS_H_
#define CURRENTFLOWGROUPCLOSENESS_H_
#include "../graph/Graph.h"
#include "Centrality.h"

namespace NetworKit {

/**
 *
 */
    class CurrentFlowGroupCloseness : public NetworKit::Centrality {
    public:
/**
 * Current Flow Group Closeness Centrality finds a Subset S of k vertices with an approximation
 * of (1-k/(k-1)*1/e-e) to the optimum in nearly linear runtime for any e>0.
 *
 * @param G Unweighted graph
 * @param k Size of the group of nodes
 * (-------------------------------------@param L Laplacian matrix of G-----------------------------------)
 * @param e Approximation error
 */

        CurrentFlowGroupCloseness(const Graph& G, count k = 1, double e=0);

        double w_min;
        double w_max;
        int n;
        int m;

        /**
	    *
	    */
        void run();
        /**
	    *
	    */
        /*vector GainsEST()*/
        std::vector<double> GainsEst(Graph& G, std::vector<int> S, double e);
        /**
	    *
	    */
        void run();
        /**
         *
         */
        void Generate_Gaussian_matrices(double G[][]);
        /*
         *
         */
        void matrix_delete_row_and_column(double G[][], int i);
        /*
         *
         */
        std::vector<double> ERSumEst(Graph& G,double L[][], double e);
        /*
         *
         */
        bool check(std::vector<int> S, int u);
    };

} /* namespace NetworKit */
#endif /* CURRENTFLOWGROUPCLOSENESS_H_ */
