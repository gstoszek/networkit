/*
 * LaplacianCentrality.h
 *
 *  Created on: 08.03.2018
 *      Author: Kolja Esders
 */

#ifndef LAPLACIANCENTRALITY_H_
#define LAPLACIANCENTRALITY_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 * Computes the Laplacian centrality of the graph.
 */
class LaplacianCentrality: public Centrality {
public:
	/**
	 * Constructs a LaplacianCentrality object for the given Graph @a G.
	 *
	 * @param[in] G The graph.
	 */
	LaplacianCentrality(const Graph& G);

	/**
	 * Computes laplacian centrality on the graph passed in constructor.
	 */
	void run() override;
};

} /* namespace NetworKit */
#endif /* LAPLACIANCENTRALITY_H_ */