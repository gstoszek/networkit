/*
 * CurrentFlowGroupCloseness.h
 *
 *  Created on:
 *      Author: GSTOSZEK
 */

#include "CurrentFlowGroupCloseness.h"
#include "Centrality.h"

/**
 * CSRMatrix.h for Laplacian
 */

namespace NetworKit {
/***SAMPLING*** (Number 2**/
 CurrentFlowGroupCloseness::CurrentFlowGroupCloseness(const Graph &G, count k, count e)
 {

    }
    /**
     *
     */
    void CurrentFlowGroupCloseness::run() {
        /*Graph to Laplacian*/
        int L[n][n];
        /*initialize ERSumsEst*/
        std::vector<double> r= ERSumEst(G, L, e / 3);
        /*Group S*/
        std::vector<int> S;
        /*find argmin of u*/
        double current_min;
        int current_vertice;
        current_min = r.at(0);
        for(int i=1; i<n; i++){
            if(current_min > r.at(i)){
                   current_vertice = i;
            }
        }
        S.push_back(current_vertice);

        /*repeat GainsEst till k*/
        std::vector<double> g = GainsEst(G,S,e);
        for(int i= 1; i<k; i++){
            /*step a*/
            g = GainsEst(G,S,e/2);
            /*step b*/
            /*!!!!!!!!!!!!!!!!*/
            current_min = g.at(0);
            for(int i=1; i<n; i++){
                if(check(S,i) = true ) {
                    if (current_min > g.at(i)) {
                        current_vertice = i;
                    }
                }
            }
            S.push_back(current_vertice);
        }
        return S;
    }
    /**
     *
     * @param G
     * @param S
     * @param e
     * @return
     */
    std::vector<double> CurrentFlowGroupCloseness::GainsEst(Graph &G, std::vector<int> S, double e)
    {
        /* Set delta_1, delta_2 and delta_3*/
        double delta_1;
        double delta_2;
        double delta_3;

        int p;
        int q;
        int r;

        /*check*/
        delta_1 = (w_max*e)/(27*w_min)*std::sqrt(((1-e/9)/(1+e/9)*n));

        /*check*/
        delta_2 = (1/4*std::pow(n,4)*w_max)*std::sqrt(e*std::sqrt(std::pow(w_min,3))/(9*n));

        delta_3 = delta_2;

        /* Generate random Gaussian matrices P,Q,R */
        /*check */
        p = std::ceil(4*(std::pow(e/9,2)/2-1/(std::pow(e/9,3)/3)*ln));
        q=p;
        r=p;

        double P[p][n];
        double Q[q][m];
        double R[r][n];
        double X[n][n];

        /* MatrixtoLaplacian */
        /* !!!!!!!!!!!!!!!!! */
        X= L;
        Generate_Gaussian_matrices(P);
        Generate_Gaussian_matrices(Q);
        Generate_Gaussian_matrices(R);

        /* Let X be Diag(L_Si 1) */
        /* !!!!!!!!!!!!!!!!!!!!! */
        for(int i= 0; i<= S.end();i++){
            matrix_delete_row_and_column(X, S.at(i));
        }



    }

    /**
     *
     * @param GK
     */
    void CurrentFlowGroupCloseness::Generate_Gaussian_matrices(double GAG[][])
    {
        /*intialising standard deviation to 1.0*/
        double sigma = 1.0;
        double r, s = 2.0 * sigma * sigma;

        /* sum is for normalization*/
        double sum = 0.0;

        /*matrix orientation*/


        for (int i = 0; x < G.row ; i++)
        {
            for(int j = 0; y < G.column; j++)
            {
                r = sqrt(i*i + j*j);
                G[i][j] = (exp(-(r*r)/s))/(M_PI * s);
            }
        }
    }
    /*
     *
     */
    void CurrentFlowGroupCloseness::matrix_delete_row_and_column(double G[][],i){

        for(int k=0; k<n; k++){
            G[k][i]= -1;
            G[i][k]= -1;
        }
    }

    std::vector<double> CurrentFlowGroupCloseness::ERSumEst(Graph &G, L[][], double e) {
        /*Set delta */
        double delta;
        delta = e/(9*std::pow(n,2))*std::sqrt(((1-e/3)*w_min)/((1+e/3)*w_max));

        /*Generate a random Gaussian matrix Q*/
        int q;
        q = std::ceil(4*(std::pow(e/3,2)/2-1/std::pow(e/3,3)/3));
        double Q[q][n];

        /*Compute Q sqrt(W) B*/
        /* !!!!!!!!!!!!!!!!!!*/

        /*Solver for approximation of Z*/

    }

    bool CurrentFlowGroupCloseness::check(std::vector<node> S, int u) {
        bool check = true;
        for(int i = 0; i< S.size();i++){
            if(S.at(i)=u){
                check= false;
            }
        }
        return check;
    }
}