/*
 * CurrentFlowGroupCloseness.h
 *
 *      Author: gstoszek
 */

#ifndef CURRENTFLOWGROUPCLOSENESS_H_
#define CURRENTFLOWGROUPCLOSENESS_H_

#include "Centrality.h"
#include "../algebraic/CSRMatrix.h"
#include "../numerics/LAMG/Lamg.h"
#include "EffectiveResistanceDistance.h"
#include "ERD2.h"
#include <armadillo>

namespace NetworKit {
    /*
     *
     */
    class CurrentFlowGroupCloseness: public NetworKit::Centrality {

    public:

        CurrentFlowGroupCloseness(const Graph& G, const double epsilon=0.1, const double delta=0.1, const double universalConstant=1.0, const int groupsize = 2);

        void run();

        void setCFGCC(double newCFGCC);
        std::vector<node> getNodesofGroup();

        double getCFGCC();


    private:

        double epsilon;
        double delta;
        count r;
        double universalConstant;
        count groupsize;
        std::vector<node> S;
        std::vector<std::vector<double>> D;

        double CFGCC;

    };


} /* namespace NetworKit */
#endif /* CURRENTFLOWGROUPCLOSENESS_H_ */
