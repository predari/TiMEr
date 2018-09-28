#include <argtable2.h>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 

#include "./../KaHIPv0.62/lib/data_structure/graph_access.h"
#include "./../KaHIPv0.62/lib/io/graph_io.h"
#include "./../KaHIPv0.62/lib/tools/macros_assertions.h"
#include "./parse_parameters.h"
#include "./../KaHIPv0.62/lib/partition/graph_partitioner.h"
#include "./../KaHIPv0.62/lib/partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "./../KaHIPv0.62/lib/tools/quality_metrics.h"
#include "./../KaHIPv0.62/lib/tools/random_functions.h"
#include "./../KaHIPv0.62/lib/tools/timer.h"

int main(int argn, char **argv) {
  PartitionConfig partition_config, repartition_config;
  std::string graph_filename;
  
  bool is_graph_weighted = false;
  bool suppress_output   = false;
  bool recursive         = false;
  
  int ret_code = parse_repart_parameters(argn, argv, 
					 partition_config, 
					 repartition_config, 
					 graph_filename, 
					 is_graph_weighted, 
					 suppress_output, recursive);
  if(ret_code) {
    return 0;
  }
  
  std::streambuf* backup = std::cout.rdbuf();
  std::ofstream ofs;
  ofs.open("/dev/null");
  if(suppress_output) {
    std::cout.rdbuf(ofs.rdbuf()); 
  }
  partition_config.LogDump(stdout);

  printf("  partition_config: k and new_k are %d and %d\n", partition_config.k, partition_config.new_k);
  printf("repartition_config: k and new_k are %d and %d\n", repartition_config.k, repartition_config.new_k);
  printf("  partition_config: imbalance and new_imbalance are %f and %f\n", partition_config.imbalance, partition_config.new_imbalance);
  printf("repartition_config: imbalance and new_imbalance are %f and %f\n", repartition_config.imbalance, repartition_config.new_imbalance);





  std::cout.rdbuf(backup);
  return 0;
}
