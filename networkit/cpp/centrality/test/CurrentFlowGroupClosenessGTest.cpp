/*
 * CurrentFlowGroupClosenessGTest.h
 *
 *      Author: gstoszek
 */


#include "CurrentFlowGroupClosenessGTest.h"
#include "../CurrentFlowGroupCloseness.h"
#include "../../generators/ErdosRenyiGenerator.h"

namespace NetworKit {

    TEST_F(CurrenFlowGroupClosenessGTest, CurrenFlowGroupClosenessGTest) {

    /* Graph:
           0
          / \
         1 - 5
        / \
       2   3 - 4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(3, 4);
    G.addEdge(0, 5);
    G.addEdge(1, 5);

    }//CurrentFlowGroupClosenessGTEST

    TEST_F(CurrenFlowGroupClosenessGTest, CurrenFlowGroupClosenessGTest2) {

    /* Graph:
           0
          / \
         1 - 5
        / \   \
        2  3 - 4
    */
    count n = 6;
    Graph G(n);

    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(3, 4);
    G.addEdge(0, 5);
    G.addEdge(1, 5);

    /**new EDGE(4,5)**/
    G.addEdge(4, 5);

    }

}