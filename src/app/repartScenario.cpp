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

#include "parse_parameters.h"
#include "quality_metrics.h"
#include "./../lib/repartitioning/graph_repartitioner.h"
#include "configuration.h"
#include "partition/graph_partitioner.h"
#include "./../lib/build_app_graph.h"

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

void repartition(graph_access & G, PartitionConfig &repartition_config) {
  graph_repartitioner repartitioner;
  quality_metrics qm;

  repartition_config.use_fullmultigrid = false; //TODO: is this ok?                                 
  repartition_config.use_wcycles = false;       //TODO: is this ok?                                 

  repartition_config.label_propagation_refinement = false;
  repartition_config.quotient_graph_refinement_disabled = true;
  repartition_config.corner_refinement_enabled = false;
  repartition_config.softrebalance = false;
  repartition_config.refinement_type = REFINEMENT_TYPE_FM;//REFINEMENT_TYPE_FM_FLOW;                

  //repartition_config.diffusion_based_repartitioning = true;
  repartition_config.diffusive_refinement = true;
  repartition_config.trunccons_limited = true;
  repartition_config.trunccons_numDiffIters = 14;
  repartition_config.trunccons_numConsols = 9;
  repartition_config.trunccons_limited_useAllNeighbors = true;
  repartition_config.trunccons_limited_sizeOfNeighborhood = 8;

  repartition_config.scale_balancing = true;
  repartition_config.scale_balancing_iters = 16;
  repartition_config.flow_based_balancing = false;
  repartition_config.greedy_balancing = true;





  timer t;
  t.restart();
  repartitioner.perform_partitioning(repartition_config, G);
  /*complete_boundary boundary(&G);
  boundary.build();
  cycle_refinement cr;
  cr.perform_refinement(partition_config, G, boundary);*/
  std::cout <<  "time spent for repartitioning " << t.elapsed()  << std::endl;
       
  // output some information about the partition that we have computed 
  std::cout << "cut \t\t"         << qm.edge_cut(G)                 << std::endl;
  std::cout << "finalobjective  " << qm.edge_cut(G)                 << std::endl;
  std::cout << "bnd \t\t"         << qm.boundary_nodes(G)           << std::endl;
  std::cout << "balance \t"       << qm.balance(G)                  << std::endl;
  std::cout << "finalbalance \t"  << qm.balance(G)                  << std::endl;
  std::cout << "max_comm_vol \t"  << qm.max_communication_volume(G) << std::endl;
}

// void fromUnbalanced(const char* file, int k, int genImbalance, int reImbalance, int mode, const char*  file_out) {
//   std::cout << "File: \t\t\t\t" << file << std::endl;
//   std::cout << "#Partitions: \t\t\t" << k << std::endl;
//   std::cout << "generator Imbalance: \t\t" << genImbalance << std::endl;
//   std::cout << "repartitioner Imbalance: \t" << reImbalance << std::endl;
//   std::cout << "mode: \t\t\t\t" << mode << std::endl;

//   PartitionConfig in_config;
//   graph_access G;
//   generateInput(in_config, G, k, genImbalance, file, mode);

//   repartition(G, k, reImbalance, mode);

//   if (file_out != NULL) graph_io::writePartition(G, file_out);
// }

int readAncestors(int* ancestors, graph_access & H, std::string filename_ancestors) {
  FILE* f = fopen(filename_ancestors.c_str(), "r");
  for (NodeID i = 0; i < H.number_of_nodes(); ++i) {
    fscanf(f, "%d", &ancestors[i]);
  }
  fclose(f);
  return(0);
}

int main(int argn, char** argv) {
  PartitionConfig partition_config, repartition_config;
  std::string graph_filename;
  graph_access H; //graph to be repartitioned
  graph_access G; //old graph (optional)
  //parsing the parameters
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
    return(0);
  }
  
  printf("  partition_config: k and new_k are %d and %d\n", partition_config.k, partition_config.new_k);
  printf("repartition_config: k and new_k are %d and %d\n", repartition_config.k, repartition_config.new_k);
  printf("  partition_config: imbalance and new_imbalance are %f and %f\n", partition_config.imbalance, partition_config.new_imbalance);
  printf("repartition_config: imbalance and new_imbalance are %f and %f\n", repartition_config.imbalance, repartition_config.new_imbalance);

  //building H
  buildAppGraph(partition_config, graph_filename, H);

  //if there is no G
  if(repartition_config.filename_old_graph == "") {
    partition(H, partition_config);
    repartition(H, repartition_config);
    return(0);  
  }

  //building G (possibly partitioned)
  buildAppGraph(partition_config, repartition_config.filename_old_graph, G);

  if(partition_config.graph_allready_partitioned  == false) {
    //partition G
    partition(G, partition_config);
  } else {
    //transplant partition of G onto H
    int* ancestors = new int[H.number_of_nodes()];
    readAncestors(ancestors, H, repartition_config.filename_ancestors);
    forall_nodes(H, i) {
      H.setPartitionIndex(i, G.getPartitionIndex(ancestors[i]));
    } endfor    
  }

  //repartition H
  repartition(H, repartition_config);

  // // write the repartition to the disc                                            
  // std::stringstream filename;
  // if(!repartition_config.filename_output.compare("")) {
  //   // no output filename given                                           
  //   filename << "tmppartition" << repartition_config.k;
  // } else {
  //   filename << repartition_config.filename_output;
  // }
  // //printf("filename is %s\n", (filename.str()).c_str());
  // graph_io::writePartition(H, filename.str());

  return 0;
}
