#include "graph_repartitioner.h"

#include "coarsening/coarsening.h"
#include "initial_partitioning/initial_partitioning.h"
#include "uncoarsening/uncoarsening.h"
#include "w_cycles/wcycle_partitioner.h"
#include "diffusive_uncoarsening.h"

//TODO: call new refinement methods

void graph_repartitioner::single_run( PartitionConfig & config, graph_access & G) {
  for( unsigned i = 1; i <= config.global_cycle_iterations; i++) {
    PRINT(std::cout <<  "vcycle " << i << " of " << config.global_cycle_iterations  << std::endl;)
     if(config.use_wcycles || config.use_fullmultigrid)  {
       wcycle_partitioner w_partitioner;
       w_partitioner.perform_partitioning(config, G);
      } else {
	coarsening coarsen;
	initial_partitioning init_part;
	//uncoarsening uncoarsen;
	diffusive_uncoarsening uncoarsen;

	graph_hierarchy hierarchy;

	coarsen.perform_coarsening(config, G, hierarchy);
	init_part.perform_initial_partitioning(config, hierarchy);
	uncoarsen.perform_uncoarsening(config, hierarchy);
      }
    config.graph_allready_partitioned = true;
    config.balance_factor             = 0;
  }
}
