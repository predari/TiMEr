#include <argtable2.h>
#include <cmath>
#include <cstring> 
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <list>
#include <regex.h>
#include <sstream>
#include <cstdlib>
#include <queue>

#include "./parse_parameters.h"
#include "quality_metrics.h"
#include "./../lib/mapping/processor_graph_access.h"
#include "./../lib/mapping/build_processor_graphs.h"
#include "./../lib/build_app_graph.h"
#include "./../lib/mapping/mapping_algorithms.h"
#include "./../lib/mapping/build_processor_graphs.h"

//#include "./../lib/mapping/mapping_algorithms.h"
#include "tools/timer.h"

/******************************************************/
/*                     partition                      */
/******************************************************/
void partition(graph_access & G, PartitionConfig& partition_config) {
  timer t;
  graph_partitioner partitioner;
  quality_metrics qm;

  t.restart();
  partitioner.perform_partitioning(partition_config, G);
  cout <<  "time spent for partitioning " << t.elapsed()  << endl;
       
  // output some information about the partition that we have computed 
  cout << "cut \t\t"         << qm.edge_cut(G)                 << endl;
  cout << "finalobjective  " << qm.edge_cut(G)                 << endl;
  cout << "bnd \t\t"         << qm.boundary_nodes(G)           << endl;
  cout << "balance \t"       << qm.balance(G)                  << endl;
  cout << "finalbalance \t"  << qm.balance(G)                  << endl;
  cout << "max_comm_vol \t"  << qm.max_communication_volume(G) << endl;
}
  
/******************************************************/
/*                  buildCommGraph                    */
/******************************************************/
void buildCommGraph (graph_access & G, PartitionConfig & partition_config, graph_access & C) {

  unsigned k = partition_config.k;//number of nodes in C
  unsigned ksq = k * k;

  //find number of edges in C and edge weights
  EdgeWeight edgeWeights[ksq];
  for(unsigned i = 0; i < ksq ; i ++) {
    edgeWeights[i] = 0;
  }
  unsigned nmbEdges = 0;
  forall_nodes(G, u) {
    PartitionID uBlock = G.getPartitionIndex(u);
    forall_out_edges(G, e, u) {
      NodeID v = G.getEdgeTarget(e);
      PartitionID vBlock = G.getPartitionIndex(v);
      if(uBlock != vBlock) {
	unsigned indexC = (unsigned) ((uBlock * k) + vBlock);
	if(edgeWeights[indexC] == 0) {
	  nmbEdges++;
	}
	(edgeWeights[indexC])++;
      }
    } endfor    
  } endfor

  //build C from k, nmbEdges, edgeWeights[]
  C.start_construction((NodeID) k, nmbEdges);
  for(unsigned i = 0; i < k; i++) {
    NodeID u = C.new_node();    
    C.setNodeWeight(u, 1);
    C.setPartitionIndex(u, 0);
    for(unsigned j = 0; j < k; j ++) {
      unsigned indexC = ((i * k) + j);
      if(edgeWeights[indexC] != 0) {
	EdgeID e = C.new_edge(i, j);
	C.setEdgeWeight(e, edgeWeights[indexC]);
      }
    }
  }
  C.finish_construction();

  printf("nmbEdges is %d\n", nmbEdges);
  forall_nodes(C, u) {
    forall_out_edges(C, e, u) {
      NodeID v = C.getEdgeTarget(e);
      bool problem = true;
      forall_out_edges(C, f, v) {
	NodeID w = C.getEdgeTarget(f);
	if(w == u) {
	  problem = false;
	}
      } endfor
      if(problem == true) {
	printf("Was soll den das? KEINE RUECKKANTE!!!\n");
      }
    } endfor
  } endfor
}

/******************************************************/
/******************************************************/
/************************* main ***********************/
/******************************************************/
/******************************************************/
int main(int argn, char** argv) {
  PartitionConfig partition_config, map_config;
  string graph_filename;
  processor_graph_access P;
  graph_access G, C;

  //parsing the parameters
  bool is_graph_weighted = false;
  bool suppress_output   = false;
  bool recursive         = false;
  
  int ret_code = parse_map_parameters(argn, argv, 
				      partition_config, 
				      map_config, 
				      graph_filename, 
				      is_graph_weighted, 
				      suppress_output, recursive);
  if(ret_code) {
    error("Parsing of parameters failed.", 600);
  }

  //building application graph G (possibly partitioned)
  timer t;
  t.restart();
  buildAppGraph(partition_config, graph_filename, G);
  if(partition_config.input_partition == "") {
    cout << "Building the application graph took " << t.elapsed() << " seconds \n";
  } else {
    cout << "Building the partitioned application graph took " << t.elapsed() << " seconds \n";
  }

  //if G is not yet partitioned
  if(partition_config.input_partition == "") {
    t.restart();
    partition(G, partition_config);
    cout << "Partitioning the application graph took " << t.elapsed() << " seconds \n";
  }

  //building communication graph C from partitioned G
  t.restart();
  buildCommGraph(G, partition_config, C);
  cout << "Building the communication graph took " << t.elapsed() << " seconds \n";

  if(partition_config.filename_output != "") {
    /* Write C to the disk */
    graph_io::writeGraphWeighted(C, partition_config.filename_output);
  }
  
  return 0;
}
