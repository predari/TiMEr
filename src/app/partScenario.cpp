#include <argtable2.h>
#include <cstdio>
#include <cmath>
#include <cstring> 
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <list>
#include <regex.h>
#include <sstream>
#include <string>

#include "./parse_parameters.h"
#include "data_structure/graph_access.h"
#include "./../lib/build_app_graph.h"
#include "quality_metrics.h"
#include "../lib/karma_config.h"
#include "configuration.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "partition/graph_partitioner.h"
#include "random_functions.h"
#include "timer.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"

void partition(graph_access & G, PartitionConfig& partition_config) {
  timer t;
  graph_partitioner partitioner;
  quality_metrics qm;

  t.restart();
  partitioner.perform_partitioning(partition_config, G);
  std::cout <<  "time spent for partitioning " << t.elapsed()  << std::endl;
       
  // output some information about the partition that we have computed 
  std::cout << "cut \t\t"         << qm.edge_cut(G)                 << std::endl;
  std::cout << "finalobjective  " << qm.edge_cut(G)                 << std::endl;
  std::cout << "bnd \t\t"         << qm.boundary_nodes(G)           << std::endl;
  std::cout << "balance \t"       << qm.balance(G)                  << std::endl;
  std::cout << "finalbalance \t"  << qm.balance(G)                  << std::endl;
  std::cout << "max_comm_vol \t"  << qm.max_communication_volume(G) << std::endl;
}

int main(int argn, char** argv) {
  PartitionConfig partition_config;
  std::string graph_filename;
  graph_access G;

  //parsing the parameters
  bool is_graph_weighted = false;
  bool suppress_output   = false;
  bool recursive         = false;
  
  int ret_code = parse_parameters(argn, argv, 
				  partition_config, 
				  graph_filename, 
				  is_graph_weighted, 
				  suppress_output, recursive);
  if(ret_code) {
    return(0);
  }
  
  printf("partition_config: k and imbalance are %d and %f\n", partition_config.k, partition_config.imbalance);

  //building G
  buildAppGraph(partition_config, graph_filename, G);
  partition(G, partition_config);

  // write the partition to the disc                                            
  std::stringstream part_filename_stream;
  std::string part_filename;
  if(!partition_config.filename_output.compare("")) {
    // no output filename given
    part_filename = graph_filename;
    part_filename.erase(part_filename.end()-6, part_filename.end());
    part_filename = part_filename + "_" + std::to_string(partition_config.k) + ".partition";
    part_filename_stream << "tmppartition" << partition_config.k;
  } else {
    part_filename_stream << partition_config.filename_output;
  }
  //graph_io::writePartition(G, part_filename);//RG version
  graph_io::writePartition(G, part_filename_stream.str());//KaHIP version

  return 0;
}
