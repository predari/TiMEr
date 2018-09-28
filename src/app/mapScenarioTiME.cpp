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
#include <unordered_map>
#include <queue>

#include "./parse_parameters.h"
#include "../lib/mapping/alma_graph_access.h"
#include "quality_metrics.h"
#include "./../lib/build_app_graph.h"
#include "./../lib/mapping/processor_graph_access.h"
#include "./../lib/mapping/build_processor_graphs.h"
#include "./../lib/mapping/mapping_algorithms.h"
#include "../KaHIPv0.62/lib/io/graph_io.h"
#include "tools/timer.h"
#include <cmath>

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
/*                      parseString                   */
/******************************************************/
int parseString(string & str, vector<int> & vec) {
  bool finished = false;
  while (!finished) {
    int strPos = str.find('x');
    if (strPos == ((int)(string::npos))) {
	finished = true;
	vec.push_back(atoi(str.substr(0, str.size()).c_str()));     	   
    } else {
      vec.push_back(atoi(str.substr(0, strPos).c_str()));
      str = str.substr(strPos+1, str.size()-strPos);
    }
    
    if(vec.back() < 1) error("No negative or zero weights.", 220);
  }
  return(vec.size());
}

/******************************************************/
/*                  buildProcGraph                    */
/******************************************************/
NodeID buildProcGraph(processor_graph_access & P, PartitionConfig & map_config) {
  vector<int> extensions;//e.g., '16x16' for 2D grid or '4x3x3' for tree, depth 3, top-down
  vector<int> bandwidths;//e.g., '3x5x7x9' for 16x16 torus or tree of depth 4, top-down
  int dimORdepthExt = parseString(map_config.proc_graph_extensions, extensions);
  int dimORdepthBand = parseString(map_config.proc_graph_bandwidths, bandwidths);

  //calculate number of compute nodes (of P)
  NodeID numberOfComputeNodes = 1;
  for(int i = 0; i < dimORdepthExt; i++) {
    numberOfComputeNodes *= (NodeID) (extensions[i]);
  }

  //consistency checks
  if(dimORdepthExt < 2) {
   error("Extension, respectively depth of processor graph must be at least two.", 200);
  }
  int size  = extensions[0];
  if((map_config.proc_graph_type != "tree") &&(map_config.proc_graph_type != "hypercube")) {//i.e., processor graph is a grid or a torus
    //Size is the extension of grid or torus in first direction.
    //For the time being, the extension should be the same in all directions.
    //The extension (size) should a power of two, exponent >= 1
    int exponent = 0;
    int sizee = size;//check whether size is a power of two
    while(sizee > 1) {
      if((sizee / 2) * 2 == sizee) {
	sizee/=2;
	exponent++;
      } else {
	error("Extension of grids and tori must a power of two.", 203);
	sizee = 0;
      }
    }
    if(dimORdepthExt > 3) {
      error("Extension of grids and tori must not be greater than three.", 201);
    }
    for(int i = 1; i < dimORdepthExt; i++) {
      if(extensions[i] != size) {
	error("Grids and tori must have same extensions in all dimensions.", 202);
      }
    }
    if(dimORdepthBand != exponent) {
      error("Vector of bandwidths should have length equal to log_2(extension) (for grids and tori).", 205);
    }
  } else if (map_config.proc_graph_type == "tree") {
    if(dimORdepthExt != dimORdepthBand) {
      error("Vector of extensions and vector of bandwidths must have same length (for trees).", 204);
    }
  } else if (map_config.proc_graph_type == "hypercube") {
    if((dimORdepthBand != 1) || (bandwidths[0] != 1)) {
      error("Hypercube bandwidth should be 1 (single number).", 501);
    }
    for(int i = 0; i < dimORdepthExt; i++) {
      if(extensions[i] != 2) {
	error("Extensions for hypercube should all be equal to 2", 503);
      }
    }
  } 
  
  //distinguish between different types of processor graphs
  //printf("map_config.proc_graph_type and dimORdepthExt are %s and %d\n", map_config.proc_graph_type.c_str(), dimORdepthExt);
  if(map_config.proc_graph_type == "grid") {
    if(dimORdepthExt == 2) {
      build2DGrid(P, size, bandwidths);
    } else {
      build3DGrid(P, size, bandwidths);
    }
  } else if (map_config.proc_graph_type == "torus") {
    if(dimORdepthExt == 2) {
      build2DTorus(P, size, bandwidths);
    } else {
      build3DTorus(P, size, bandwidths);
    }
  } else if (map_config.proc_graph_type == "tree") {
    buildFatTree(P, extensions, bandwidths);
  } else if (map_config.proc_graph_type == "hypercube") {
    buildHypercube(P, dimORdepthExt);
  } else {
    error("This type of processor graph does not exist.", 206);
  }
  return(numberOfComputeNodes);
}

/******************************************************/
/*                  buildCommGraph                    */
/******************************************************/
void buildCommGraph(graph_access & G, PartitionConfig & partition_config, graph_access & C) {

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
      //cout << "uBlock and vBlock are " << uBlock << " and " << vBlock << endl;
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
}

/******************************************************/
/*                   checkMapping                     */
/******************************************************/
bool checkMapping(processor_graph_access & P, PartitionID k, vector<int> & mapping) {
  // for(PartitionID i = 0; i < k; i++) {
  //   printf("i and mapping[i] is %d and %d\n",i, mapping[i]);
  // }
  bool valid = true;
  NodeID kp = P.number_of_nodes();
  vector<int> numbHits(kp, 0);
  for(PartitionID i = 0; i < k; i++) {
    (numbHits[mapping[i]])++;
  }

  //Mapping is correct if and only if
  //
  //(i)   every compute node is mapped onto,
  //(ii)  none of the P's nodes is mapped onto more than once, and
  //(iii) none of the switches is mapped onto.
  //
  //Note that (i) implies (ii), (iii),
  //and that the weight of a node of P is 1 [0] if
  //the node is a compute node [switch]. 
  for(NodeID j = 0; j < kp; j++) {
    if(numbHits[j] != 1) {
      if((numbHits[j] == 0) && (P.getNodeWeight(j) == 1)) {
	error("Mapping not surjective.", 130);
	valid = false;
      }     
      if(numbHits[j] > 1)
	error("Mapping not injective.", 131);
      valid = false;
    } else {//i.e., numbHits[j] == 1
      if (P.getNodeWeight(j) != 1) {
	error("Process mapped to switch.", 132);	
      }
    }
  }
  return(valid);
}

/******************************************************/
/*                  evaluateMapping                   */
/******************************************************/
int evaluateMapping(graph_access & C, processor_graph_access & P, double mappingDuration, vector<int> & mapping, unsigned int call) {

  //check if mapping is valid
  PartitionID k = C.number_of_nodes();
  bool valid = checkMapping(P, k, mapping);
  if(valid == false) {
    return(0);
  }

  //compute maximum congestion, maximum dilation, average dilation
  int sumDilation = 0;
  int maxDilation = 0;
  vector<int> congestion(P.number_of_edges(), 0);
  int maxCongestion = 0;


  // checking that entries in distanceMatrix does not exceed the diameter of the processors graphs
  // maria
  // BECAREFULL hard coded for the 3d grid 8*8*8 case!
  // int Pnodes = P.number_of_nodes();
  // int DiameterGrid3D = 21; // 3*(std::pow (Pnodes,1/3) -1);
  // for (int k = 0; k < Pnodes; k++) {
  //   for (int i = 0; i < Pnodes; i++) {
  //     if( P.distanceMatrix[k][i] > DiameterGrid3D )
  // 	printf("Shortest path between %d and %d is %d > %d!\n",k,i,P.distanceMatrix[k][i], DiameterGrid3D);      
  //   }
  // }
  // maria
  for(PartitionID i = 0; i < k; i++) {	
    forall_out_edges(C, edgeC, i) {
      if(mapping[i] <= mapping[C.getEdgeTarget(edgeC)]) {//only one edge direction considered
	//find dilation of edgeC, update sumDilation, maxDilation
	unsigned int start = mapping[i];
	unsigned int target = mapping[C.getEdgeTarget(edgeC)];
	int currDilation = (P.distanceMatrix[start][target]) * (C.getEdgeWeight(edgeC));
	sumDilation += currDilation;
	if(currDilation > maxDilation) {
	  maxDilation = currDilation;
	}
	
	//update congestion[]
	unsigned int current = target;
	unsigned int next = target;
	while (current != start) {
	  current = P.predecessorMatrix[start][current];
	  if (next >= current) {
	    forall_out_edges(P, edgeP, current) {
	      if(P.getEdgeTarget(edgeP) == next) {
		congestion[edgeP] += C.getEdgeWeight(edgeC);
	      }
	    } endfor
          } else {
	    forall_out_edges(P, edgeP, next) {
	      if (P.getEdgeTarget(edgeP) == current) {
		congestion[edgeP] += C.getEdgeWeight(edgeC);
	      }
	    } endfor
	  }
	  next = current;
	}
      }
    } endfor
  } 

  forall_edges(P, edgeP) {
    //printf("edge %d: traffic=%d , weight=%d\n",edgeP, congestion[edgeP], P.getEdgeWeight(edgeP));
    (congestion[edgeP]) /= P.getEdgeWeight(edgeP);//recall that edge weight indicates bandwidth
  } endfor

  forall_edges(P, edgeP) {
    if (congestion[edgeP] > maxCongestion) {
      maxCongestion = congestion[edgeP];
      printf("The edge with the maxCongestion is %d\n",edgeP);
    }
  } endfor

  double avgDilation = ((double) sumDilation) / ((double) (C.number_of_edges() / 2));

  cout << call << ":::Mapping duration: " << mappingDuration << "\n";
  cout << call << ":::Maximum congestion: " << maxCongestion << "\n";
  cout << call << ":::Maximum dilation: " << maxDilation << "\n";
  cout << call << ":::Average dilation: " << avgDilation << "\n";
  // maria cout << call << ":::Communication Graph: meanDegree: " << C.getMeanDegree(C) << "minDegree:" << C.getMinDegree(C) << "maxDegree:" << C.getMaxDegree(C) << "\n";
  
  
  return(0);
}

/******************************************************/
/******************************************************/
/************************* main ***********************/
/******************************************************/
/******************************************************/
int main(int argn, char** argv) {
  PartitionConfig partition_config, map_config;
  string graph_filename;

  alma_graph_access g;
  processor_graph_access P;
  graph_access C;

  //parsing the parameters
  bool is_graph_weighted = false;
  bool suppress_output   = false;
  bool recursive         = false;

  printf(" Input \'partition\' before process \n");
  for(unsigned int i = 0; i < 15; i++) {
    printf("%d ", partition_config.input_partition[i]);
  }
  
  int ret_code = parse_map_parameters(argn, argv, 
				      partition_config, 
				      map_config, 
				      graph_filename,
				      is_graph_weighted, 
				      suppress_output, recursive);
  if(ret_code) {
    error("Parsing of parameters failed.", 600);
  }
  printf(" Input \'partition\' after parsing \n");
  for(unsigned int i = 0; i < 15; i++) {
    printf("%d ", partition_config.input_partition[i]);
  }
  
  printf("map_config: mapping_algo and type of processor graph are %s and %s\n",
	 map_config.mapping_algo.c_str(), map_config.proc_graph_type.c_str());



  //building processor graph P
  timer t;
  NodeID numberOfComputeNodes = buildProcGraph(P, map_config);
  cout << "Building the processor graph took " << t.elapsed() << " seconds \n";
  //graph_io::writeGraphWeighted(P, "mesh_16x16.graph");
  
  //building application graph g (possibly partitioned)
  t.restart();
  buildAppGraph(partition_config, graph_filename, g);
  if(partition_config.input_partition == "") {
    cout << "Building the application graph took " << t.elapsed() << " seconds \n";
  } else {
    cout << "Building the partitioned application graph took " << t.elapsed() << " seconds \n";
  }
  cout << "Application Graph --> meanDegree: " << g.getMeanDegree(g) << "   minDegree:" << g.getMinDegree(g) << "   maxDegree:" << g.getMaxDegree(g) << "\n";

  //if g is not yet partitioned
  if(partition_config.input_partition == "") {
    cout << "======================== Partition ===========================" << endl;
    t.restart();
    partition(g, partition_config);
    cout << "Partitioning the application graph took " << t.elapsed() << " seconds \n";
  }

  // Evaluate the initial partition
  quality_metrics qm;       
  // output some information about the partition that we have computed
  cout << "evaluating the initial partition" << endl;
  cout << "cut \t\t"         << qm.edge_cut(g)                 << endl;
  cout << "finalobjective  " << qm.edge_cut(g)                 << endl;
  cout << "bnd \t\t"         << qm.boundary_nodes(g)           << endl;
  cout << "balance \t"       << qm.balance(g)                  << endl;
  cout << "finalbalance \t"  << qm.balance(g)                  << endl;
  cout << "max_comm_vol \t"  << qm.max_communication_volume(g) << endl;

  //building communication graph C from partitioned g
  t.restart();
  buildCommGraph(g, partition_config, C);
  cout << "Building the communication graph took " << t.elapsed() << " seconds \n";
  if(numberOfComputeNodes != C.number_of_nodes()) {
    error("The processor graph and communication graph have different numbers of nodes", 400);
  }
  //mapping
  cout << "======================== MAPPING ("<< map_config.mapping_algo << ") ===========================" << endl;

  t.restart();
  vector<int> mapping;//to encode mapping C --> P
  // temporarily comment the following (4) lines
  bool OK = map(g, C, P, map_config, mapping);
  if(OK == false) {
    return 0;
  }
  // the mapping coresponds to the initial partition as given from the input partition 
  // for(unsigned int i = 0; i < C.number_of_nodes(); i++) {
  //   mapping.push_back( partition_config.input_partition[i]);
  // }

  printf("Initial \'partition\'\n");
  for(unsigned int i = 0; i < 15; i++) {
    printf("%d ", partition_config.input_partition[i]);
  }
  printf("\n");
  printf("Initial mapping (should be the same as the partition)\n");
  for(unsigned int i = 0; i < 15; i++) {
    printf("%d ", mapping[i]);
  }
  printf("\n");

  double mappingDuration = t.elapsed();
  cout << "DONE WITH MAPPING" << endl;
  cout << "======================== EVALUATION MAPPING ===========================" << endl;

  //(first) evaluation of mapping
  t.restart();
  evaluateMapping(C, P, mappingDuration, mapping, 1);
  cout << "First evaluation of mapping took " << t.elapsed() << " seconds \n" << endl;

  /*******************************************************************************************/
  /******************************* POSTPROCESSING WITH ALMA **********************************/
  /*******************************************************************************************/
  cout << "======================== ALMA ===========================" << endl;
  //adapt g to new mapping
  forall_nodes(g, v) {
    g.setPartitionIndex(v, mapping[g.getPartitionIndex(v)]);    
  } endfor

  if(map_config.numberHierarchies > 0) {
    t.restart();
    int64_t initialCut = UNDEFINED_NODE;
    int64_t finalCut = UNDEFINED_NODE;
    int64_t initialTCV = UNDEFINED_NODE;
    int64_t finalTCV = UNDEFINED_NODE;
    OK = almaMapping(g, P, map_config, mapping,
		     initialCut, finalCut, initialTCV, finalTCV);
    if(OK == false) {
      return 0;
    }
    double mappingDuration = t.elapsed();
    
    t.restart();
    //First derive new partition (block IDs) of g from Hamming labels of g.
    //At the same time, determine the bijection between block IDs and block labels in g
    //Then, get the bijection between block IDs and block labels in P
    //Finally, build communication graph (a posteriori) and write the mapping
    
    //create and write the vector blockID2blockLabel
    vector<int64_t> blockID2blockLabel(partition_config.k);
    unordered_map<int64_t, int64_t> blockLabel2blockID;//auxiliary vector
    int64_t blockID = 0;
    forall_nodes(g, v) {
      //cout << "label of " << v << " is " << g.getLabel(v) << endl; 
      int64_t blockLabel_of_v = (g.getLabel(v) >> g.getNumDigitsForID());
      //cout << "blockLabel_of " << v << " is " << blockLabel_of_v << endl; 
      //cout << endl;
      if(blockLabel2blockID.find(blockLabel_of_v) == blockLabel2blockID.end()) {
	blockLabel2blockID[blockLabel_of_v] = blockID;
	blockID2blockLabel[blockID] = blockLabel_of_v;
	g.setPartitionIndex(v, blockID);
	blockID++;
      } else {
	g.setPartitionIndex(v, blockLabel2blockID[blockLabel_of_v]);
      }
    } endfor
	
    //cout << "blockID is " << blockID << endl;
	
    cout << "numDigitsL is " << g.getNumDigitsL() << ", numDigitsBlockID is " <<
	g.getNumDigitsBlockID() << " and numDigitsForID is " << g.getNumDigitsForID() << endl;
  
    //get the bijection between block IDs and block labels in P
    //create and write the vector blockLabel2vertex_in_P
    cout << "still alive 00" << endl; 
    unordered_map<int64_t, int64_t> blockLabel2vertex_in_P;
    forall_nodes(P, v) {
      blockLabel2vertex_in_P[(P.getHammingLabel(v))] = v; 
    } endfor
      
    cout << "still alive 11" << endl; 
    //build communication graph C from partitioned g WITH modifications from Alma
    t.restart();
    buildCommGraph(g, partition_config, C);
    cout << "Building the communication graph took " << t.elapsed() << " seconds \n";
    if(numberOfComputeNodes != C.number_of_nodes()) {
      error("Processor graph and the communication graph have different numbers of nodes", 400);
    }

    //write second mapping, i.e., mapping from Alma
    NodeID k = C.number_of_nodes();
    mapping.resize(k, 0);
    unsigned int offset = P.number_of_nodes() - k;
    for(unsigned int commNode = 0; commNode < C.number_of_nodes(); commNode++) {
      mapping[commNode] = blockLabel2vertex_in_P[blockID2blockLabel[commNode]] + offset;
      //cout << "mapping[" << commNode << "] is " << blockLabel2vertex_in_P[blockID2blockLabel[commNode]] + offset << endl;
      //cout << "mapping[" << commNode << "] is " << blockID2blockLabel[commNode] + offset << endl;
    }
  
    //evaluate mapping for the second time (after postprocessing by Alma)
    cout << "======================== EVALUATION ALMA ===========================" << endl;
    evaluateMapping(C, P, mappingDuration, mapping, 2);
    cout << "1:::Cut is " << initialCut << endl;
    cout << "2:::Cut is " << finalCut << endl;
    cout << "1:::TCV is " << initialTCV << endl;
    cout << "2:::TCV is " << finalTCV << endl;
    cout << "Second evaluation of mapping took " << t.elapsed() << " seconds \n" <<endl;
  }
  
  /******************************* Output *****************************/
  
  if(partition_config.filename_output != "") {
    /* Write mapping as function from vertex set of application graph g to node of processor graph */
    //first MANIPULATE partitioning, then write file with "partition"
    forall_nodes(g, v) {
      g.setPartitionIndex(v, mapping[g.getPartitionIndex(v)]);
    } endfor
    graph_io::writePartition(g, partition_config.filename_output);
  }
  
  return 0;
}
