/*
 * CurrentFlowGroupClosenessGTest.h
 *
 *      Author: gstoszek
 */


#ifndef CURRENFLOWGROUPCLOSENESSGTEST_H_
#define CURRENFLOWGROUPCLOSENESSGTEST_H_

#include "CurrentFlowGroupClosenessGTest.h"
#include "../CurrentFlowGroupCloseness.h"
#include "../../generators/ErdosRenyiGenerator.h"

namespace NetworKit {

    class CurrentFlowGroupClosenessGTest : public testing::Test {
    public:
        CurrentFlowGroupClosenessGTest() = default;
        virtual ~CurrentFlowGroupClosenessGTest() = default;
    };

}

#endif /* CURRENFLOWGROUPCLOSENESSGTEST_H_ */