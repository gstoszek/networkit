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

namespace NetworKit {
    /*
     *
     */
    class CurrentFlowGroupCloseness: public NetworKit::Centrality {

    public:
        /*
         *
         */
    CurrentFlowGroupCloseness(const Graph& G, const double epsilon=0.01, const double delta=0.1, const double universalConstant=1.0, const int groupsize = 10);
        /*
         *
         */
    void run() override;
    /*
     *
     */
    count numberOfSamples();
    /*
     *
     */

    private:

        double epsilon;
        double delta;
        count r;
        double universalConstant;
        int groupsize;
        std::vector<count> d;
        std::vector<std::vector<node>> L;
};

} /* namespace NetworKit */
#endif /* CURRENTFLOWGROUPCLOSENESS_H_ */
