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
#include "./../lib/tools/extended_quality_metrics.h"
#include "../lib/repartitioning/graph_repartitioner.h"
#include "./../lib/karma_config.h"
#include "configuration.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "random_functions.h"
#include "timer.h"

/********************************************************************/
/******************************* main *******************************/
/********************************************************************/
int main(int argn, char** argv) {
  PartitionConfig partition_config;
  std::string graph_filename;
  graph_access G;
  extended_quality_metrics eqm;

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
  
  //read G and partition of G
  std::cout <<  "Reading graph ..." << std::endl;
  graph_io::readGraphWeighted(G, graph_filename);                                                         
  std::cout <<  "Reading partition ..." << std::endl;
  graph_io::readPartition(G, partition_config.input_partition);
  partition_config.k = G.get_partition_count();
  std::cout <<  "Graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() / 2 <<  " edges."  << std::endl;
  std::cout <<  "Number of blocks is " <<  partition_config.k << "." << std::endl;

  //write statistics
  std::vector<double> block_balance(G.get_partition_count(),0.0);
  std::vector<NodeID> block_nodes(G.get_partition_count(),0);
  std::vector<NodeID> block_bnd(G.get_partition_count(),0);
  std::vector<EdgeWeight> block_internal_edges(G.get_partition_count(),0);
  std::vector<EdgeWeight> block_cut(G.get_partition_count(),0);
  std::vector<EdgeWeight> block_volume(G.get_partition_count(),0);
  double balance=eqm.balances(G, block_balance);
  int bnd=eqm.boundaryNodes(G, block_nodes, block_bnd);
  EdgeWeight cut=eqm.edgeCuts(G, block_internal_edges, block_cut);
  EdgeWeight mcv=eqm.communicationVolumes(G, block_volume);
  std::cout << "totalBalance \t\t\t" << balance << std::endl;
  std::cout << "totalNumberInternalNodes \t" << G.number_of_nodes() - bnd << std::endl;
  std::cout << "totalNumberBoundaryNodes \t" << bnd << std::endl;
  std::cout << "totalNumberInternalEdges \t" << (G.number_of_edges())/2 - cut << std::endl;
  std::cout << "totalEdgeCut \t\t\t" << cut << std::endl;
  std::cout << "maxCommVol \t\t\t" << mcv << std::endl;

  for(PartitionID i=0; i < G.get_partition_count(); i++) {
    std::cout << "balance of block " << i << " is \t\t\t" << block_balance[i]  << std::endl;
    std::cout << "Number of internal nodes in block " << i << " is \t" << block_nodes[i] - block_bnd[i]  << std::endl;
    std::cout << "Number of boundary nodes in block " << i << " is \t" << block_bnd[i] << std::endl;
    std::cout << "Number of internal edges in block " << i << " is \t" << block_internal_edges[i] << std::endl;
    std::cout << "Edge cut of block " << i << " is \t\t\t" << block_cut[i] << std::endl;
    std::cout << "Communication volume of block " << i << " is \t\t" << block_volume[i] << std::endl;
  }
  return 0;
}
