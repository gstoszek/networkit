/*
 * CurrentFlowGroupClosenessGTest.h
 *
 *      Author: gstoszek
 */


#include "CurrentFlowGroupClosenessGTest.h"
#include "../CurrentFlowGroupCloseness.h"
#include "../Betweenness.h"
#include "../../generators/ErdosRenyiGenerator.h"

namespace NetworKit {

    TEST_F(ApproxBetweennessGTest, benchApproxDiameterErdos) {
    Graph G1 = gen.generate();
    CurrentFlowGroupCloseness CurrentFlowGroupCloseness(G1, 0.05, 0.1,10);
    approx.run();
}

}