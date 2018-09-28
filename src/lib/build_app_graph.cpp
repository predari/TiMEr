#include "./build_app_graph.h"

/******************************************************/
/*                  buildAppGraph                     */
/******************************************************/
int buildAppGraph(PartitionConfig & partition_config, std::string & graph_filename, graph_access & G) {

  timer t;
  graph_io::readGraphWeighted(G, graph_filename);

  std::cout << "io time: " << t.elapsed()  << std::endl;

  if(partition_config.input_partition != "") {
    std::cout <<  "reading input partition" << std::endl;
    graph_io::readPartition(G, partition_config.input_partition);
    partition_config.graph_allready_partitioned  = true;
    partition_config.no_new_initial_partitioning = true;
    partition_config.k = G.get_partition_count();
  } else {
    G.set_partition_count(partition_config.k); 
    NodeWeight largest_graph_weight = 0;
    forall_nodes(G, node) {
      largest_graph_weight += G.getNodeWeight(node);
    }  
    endfor
     
      double epsilon = (partition_config.imbalance)/100.0;
    if(partition_config.imbalance == 0) {
      partition_config.upper_bound_partition      = (1+epsilon+0.01)*ceil(largest_graph_weight/(double)partition_config.k);
      partition_config.kaffpa_perfectly_balance   = true;
    } else {
      partition_config.upper_bound_partition      = (1+epsilon)*ceil(largest_graph_weight/(double)partition_config.k);
    }
    partition_config.largest_graph_weight       = largest_graph_weight;
    partition_config.graph_allready_partitioned = false;
    partition_config.kway_adaptive_limits_beta  = log(largest_graph_weight);
    
    std::cout <<  "block weight upper bound " <<  partition_config.upper_bound_partition  << std::endl;
    
  }

  srand(partition_config.seed);
  random_functions::setSeed(partition_config.seed);

  std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
  return(0);
}
