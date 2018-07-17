#include <stdio.h>
#include <cstdlib>
#include <queue>
#include <algorithm>

#include "./mapping_algorithms.h"
#include "./alma/hierarchy/hierarchy.h"
#include "tools/graph_extractor.h"
#include "random_functions.h"

/******************************************************/
/*                     refinement                     */
/******************************************************/
void refinement(graph_access & C, processor_graph_access & P, vector<int> & mapping, int numOfGroups, int degree, vector<int> & internalPosition) {
  
  for (int i = 0; i < numOfGroups; i++) {
    graph_access tempgraph;
    std::vector<int> nodes;
    std::vector<int> nodeWgt;
    std::vector<int> edges;
    std::vector<int> edgeWgt;
    std::vector<int> newName(numOfGroups * degree);
    std::vector<int> newNameOtherWay;
    int ctr = 0;
    nodes.push_back(0);
    nodeWgt.push_back(1);
    for(int j = 0; j < numOfGroups * degree; j++) {
      if(mapping[j] == i) {
	forall_out_edges(C,edge,j) {
	  if(mapping[C.getEdgeTarget(edge)] == i) {
	    edges.push_back(C.getEdgeTarget(edge));
	    edgeWgt.push_back(C.getEdgeWeight(edge));
	    ctr++;
	  }
	} endfor

	    newNameOtherWay.push_back(j);
	newName[j] = nodes.size() - 1;
	nodes.push_back(ctr);
	nodeWgt.push_back(1);
      }
    }

    for (unsigned int j = 0; j < edges.size(); j++) {
      edges[j] = newName[edges[j]];
    }

    tempgraph.build_from_metis_weighted(degree, &nodes[0], &edges[0], &nodeWgt[0], &edgeWgt[0]);
        
    int dest = 0;
    int target = 0;
    int weight = 0;

    forall_nodes(tempgraph, node) {
      forall_out_edges(tempgraph, edge, node) {
	if(tempgraph.getEdgeWeight(edge) > weight) {
	  weight = tempgraph.getEdgeWeight(edge);
	  dest = node;
	  target = tempgraph.getEdgeTarget(edge);
	}
      } endfor
	  } endfor

	      std::vector<int> internalMap;
    internalMap.push_back(dest);
    internalMap.push_back(target);
    for (int j = 0; j < degree; j++) {
      if(j != dest && j !=  target) {
	internalMap.push_back(j);
      }
    }
    int totalWgt = 1000000000;
    std::vector<int> bestPerm(degree);

    do {
      int mapWgt = 0;
      forall_nodes(tempgraph, node) {
	forall_out_edges(tempgraph, edge, node) {
	  mapWgt += tempgraph.getEdgeWeight(edge) * P.levelDist[internalMap[node]][internalMap[tempgraph.getEdgeTarget(edge)]];
	} endfor
	    } endfor
	  
		if(mapWgt < totalWgt) {
		  totalWgt = mapWgt;
		  std::copy (internalMap.begin(), internalMap.end(), bestPerm.begin());
		}
  	  
    } while ( std::next_permutation(internalMap.begin() + 2, internalMap.end()) );

    for (unsigned int k = 0; k < bestPerm.size(); k++) {
      internalPosition.at(newNameOtherWay[k]) = bestPerm[k];
    }
  }
}

/******************************************************/
/*                     coarseComm                     */
/******************************************************/
graph_access coarseComm (graph_access & G, std::vector<int> & map, int parts) {
  
  vector<int> nodes;
  vector<int> edges;
  vector<int> weights;

  //Combine graph G and partition.
  nodes.push_back(0);
  vector<int> iEdge;
  int lower = 0;
  int upper = 0;
  int counter = 0;

  for (int i = 0; i < parts; i++) {

    for (int j = 0; j < parts; j++) {
      iEdge.push_back(0);
    }
    
    for (unsigned int j = 0; j < map.size(); j++) {
      if(map[j] == i) {
	//	cout << j << " wird auf " << i << " gemappt!\n";
	lower = G.get_first_edge(j);
	upper = G.get_first_invalid_edge(j);
	for (int k = lower; k < upper; k++) {
	  iEdge[map[G.getEdgeTarget(k)]] += G.getEdgeWeight(k);
	}
      }
    }
    //Inserting iEdge in actual output arrays.
    for (int m = 0; m < parts; m++) {
      if(iEdge[m] > 0 && i != m) {
	edges.push_back(m);
	weights.push_back(iEdge[m]);
	//	cout << "Es gibt " << weights[counter] << " Kanten von " << i << " nach " << m << ", Counter ist " << counter << "\n";
	counter++;
      }
    }
    nodes.push_back(counter);
    iEdge.clear();
  }
  graph_access gc;	
  gc.build_from_metis(parts,&nodes[0],&edges[0]);

  forall_edges(gc, e) {
    gc.setEdgeWeight(e, weights[e]);
  } endfor 

      return gc;  
}

/******************************************************/
/*                   minimumNode                      */
/******************************************************/
int minimumNode(std::vector<int>* nodeAttribs) {
  int minNode = -1;
  int minValue = INT_MAX;
  for (unsigned int i = 0; i < nodeAttribs->size(); i++) {
    if((*nodeAttribs)[i] < minValue) {
      minNode = i;
      minValue = (*nodeAttribs)[i];
    }
  }
  return minNode;
}

/******************************************************/
/*                   maximumNode                      */
/******************************************************/
int maximumNode(std::vector<int>* nodeAttribs) {
  int maxNode = -1;
  int maxValue = -1;
  for (unsigned int i = 0; i < nodeAttribs->size(); i++) {
    if((*nodeAttribs)[i] > maxValue) {
      maxNode = i;
      maxValue = (*nodeAttribs)[i];
    }
  }
  return maxNode;
}

/******************************************************/
/*                       rand01                       */
/******************************************************/
double rand01(void){
  return(((double) rand()) / ((double)(RAND_MAX + 1.0)));
}

/******************************************************/
/*                        map                         */
/******************************************************/
bool map(alma_graph_access & g, graph_access & C, processor_graph_access & P, PartitionConfig & map_config,
	 int64_t numberHierarchies, std::vector<int> & mapping) {

  bool OK = true;
  std::string mappingAlgo = map_config.mapping_algo;
  if(mappingAlgo == "initial") {
    initialMapping(C, P, mapping);
  } else {
    if(mappingAlgo == "random") {
      randomMapping(C, P, mapping);
    } else {
      if(mappingAlgo == "drb") {
	drbMapping(C, P, mapping);
      } else {
	if(mappingAlgo == "rcm") {
	  rcmMapping(C, P, mapping);
	} else {
	  if(mappingAlgo == "multisimple") {
	    multiSimpleMapping(C, P, mapping);
	  } else {
	    if(mappingAlgo == "multikahip") {
	      multiKahipMapping(C, P, map_config.multikahip_preconfiguration, mapping);
	    } else {
	      if(mappingAlgo == "brandfass") {
		brandfassMapping(C, P, mapping);
	      } else {
		if(mappingAlgo == "brandfassL") {
		  brandfassLMapping(C, P, mapping);
		} else {
		  if(mappingAlgo == "hoefler") {
		    hoeflerMapping(C, P, mapping);
		  } else {
		    if(mappingAlgo == "hoeflerL") {
		      hoeflerLMapping(C, P, mapping);
		    } else {
		      if(mappingAlgo == "alma") {
			OK = almaMapping(g, C, P, map_config, numberHierarchies, mapping);
		      } else {
			error("Mapping name not found", 121);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return(OK);
}

/******************************************************/
/*                     bfsSearch                      */
/******************************************************/
int bfSearch(graph_access & G, NodeID startNode, NodeID & depth,  NodeID & lastNode) {
  unsigned int numberNodesFound = 1;
  vector<bool> found(G.number_of_nodes(), false);
  NodeID nodeAtEndOfLevel = startNode; 
  queue<NodeID> nodeQueue;

  // initialize queue
  nodeQueue.push(startNode);
  found[startNode] = true;

  //build queue and keep track of depth of the BFS
  depth = 0;
  while(nodeQueue.empty() == false) {
    NodeID i = nodeQueue.front();
    nodeQueue.pop();    
    forall_out_edges(G, edge, i) {
      NodeID target = G.getEdgeTarget(edge);
      if(!found[target]) {
	nodeQueue.push(target);
	found[target] = true;
	numberNodesFound++;
      } 
    } endfor

	if(i == nodeAtEndOfLevel) {
	  depth++;
	  nodeAtEndOfLevel = nodeQueue.back();
	}
  }

  //find node encountered last
  lastNode = nodeQueue.back();

  return(numberNodesFound);
}

/******************************************************/
/*                      bfsOrder                      */
/******************************************************/
void bfsOrder(graph_access & G, std::vector<NodeID> & nodeOrder, NodeID startNode) {
  vector<bool> found(G.number_of_nodes(), false);
  NodeID order = 0;
  nodeOrder[0] = startNode;
  queue<NodeID> nodeQueue;

  // initialize queue
  nodeQueue.push(startNode);
  found[startNode] = true;

  //build queue and keep track of depth of the BFS
  while(nodeQueue.empty() == false) {
    int i = nodeQueue.front();
    nodeQueue.pop();    
    forall_out_edges(G, edge, i) {
      NodeID target = G.getEdgeTarget(edge);
      if(!found[target]) {
	order++;
	nodeOrder[order] = target;
	nodeQueue.push(target);
	found[target] = true;
      } 
    } endfor
	}
}

/******************************************************/
/*                     rcmMapping                     */
/******************************************************/
void rcmMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping) {
  std::vector<NodeID> nodeOrderProc(P.number_of_nodes());
  std::vector<NodeID> nodeOrderComm(C.number_of_nodes());

  NodeID startNodeP = 0;
  NodeID startNodeC = 0;

  //find startnode in communication graph (startNodeC)
  //start_node is locally most eccentric
  NodeID maxDist;
  NodeID lastNode;
  NodeID numberNodesEncountered = bfSearch(C, startNodeC, maxDist, lastNode);
  if(numberNodesEncountered != C.number_of_nodes()) {
    error("Communication graph must be connected for rcmMapping.", 500);
  }
  bool goOn = true;
  while (goOn) {
    NodeID newDist;
    NodeID newLastNode;
    bfSearch(C, lastNode, newDist, newLastNode);
    if(newDist <= maxDist) {
      goOn = false;
    } else {
      maxDist = newDist;
      lastNode = newLastNode;
    }
  }
  startNodeC = lastNode;

  //find startnode in processor-graph, not too useful but takes nearly no time because distance already exists
  int nextNode = maximumNode(&((P.distanceMatrix)[startNodeP]));
  while(maximumNode(&((P.distanceMatrix)[nextNode])) > maximumNode(&((P.distanceMatrix)[startNodeP]))) {
    startNodeP = nextNode;
    nextNode = maximumNode(&((P.distanceMatrix)[nextNode]));
  }

  //order nodes in C and P through BFS starting at
  //startNodeC and  startNodeP, respectively
  bfsOrder(C, nodeOrderComm, startNodeC);
  bfsOrder(P, nodeOrderProc, startNodeP);

  //build mapping from nodeOrderComm, nodeOrderproc
  NodeID k = C.number_of_nodes();
  mapping.resize(k, 0);
  int j = 0;
  for (unsigned int i = 0; i < k; i++) {
    bool goOn = true;
    while (goOn) {
      NodeID procNode = nodeOrderProc[j];
      if(P.getNodeWeight(procNode) != 0) { 
	mapping[nodeOrderComm[i]] = procNode;
	j++;
	goOn = false;
      } else {
	j++;
      }	
    }
  }
}

/******************************************************/
/*                  initialMapping                    */
/******************************************************/
void initialMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping) {
  unsigned int offset = P.number_of_nodes() - C.number_of_nodes();
	
  for (unsigned int i = 0; i < C.number_of_nodes(); i++) {
    mapping.push_back(i + offset);
  }
}

/******************************************************/
/*                recursiveBisect                     */
/******************************************************/
void recursiveBisect(graph_access & G, std::vector<int> & vector, unsigned int start, unsigned int end) {
  graph_partitioner part;
  PartitionConfig partition_config;
  configuration config;
  graph_extractor extractor;

  config.fast(partition_config);
  partition_config.imbalance = 0;
  partition_config.epsilon = 0;    
  partition_config.graph_allready_partitioned = false;
  partition_config.k = 2;  
  partition_config.largest_graph_weight = end - start;
  partition_config.upper_bound_partition = ((end - start + 1) / 2);
  G.set_partition_count(2);
  part.perform_partitioning(partition_config, G);

  std::vector<unsigned int> left;
  std::vector<unsigned int> right;
  graph_access leftG;
  graph_access rightG;
  unsigned int weightLeft;
  unsigned int weightRight;
  extractor.extract_two_blocks(G, leftG, rightG, left, right, weightLeft, weightRight);
  
  std::vector<unsigned int> vecCopy(vector.size());
  std::copy (vector.begin(), vector.end(), vecCopy.begin());

  for (unsigned int i = start; i < end; i++) {
    if(i - start < weightLeft) {
      vector[i] = vecCopy[start + left[i - start]];
    } else {
      vector[i] = vecCopy[start + right[i - start - left.size()]];
    }
  }

  int middle = start + weightLeft; 

  if(middle - start >= 3) {
    recursiveBisect(leftG, vector, start, middle);
  } 

  if(end - middle >= 3) {
    recursiveBisect(rightG, vector, middle, end);
  }
}

/******************************************************/
/*                    drbMapping                      */
/******************************************************/
void drbMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping) {

  //recursive bisection on copy of C
  std::vector<int> commMap; 
  graph_access commCopy;
  C.copy(commCopy);
  for (unsigned int i = 0; i < C.number_of_nodes(); i++) {
    commMap.push_back(i);
    commCopy.setNodeWeight(i, 1);
  }
  recursiveBisect(commCopy, commMap, 0, C.number_of_nodes());

  //recursive bisection on copy of P
  std::vector<int> procMap;  
  graph_access procCopy;
  P.copy(procCopy);
  for (unsigned int i = 0; i < P.number_of_nodes(); i++) {
    procMap.push_back(i);
    procCopy.setNodeWeight(i, 1);
  }
  recursiveBisect(procCopy, procMap, 0, P.number_of_nodes());

  //find mapping based on commMap and procMap
  mapping.resize(C.number_of_nodes(), 0);
  int counter = 0;
  for (unsigned int i = 0; i < procMap.size(); i++) {
    if(P.getNodeWeight(procMap[i])) {
      mapping[commMap[counter]] = procMap[i];
      counter++;
    }
  }
}

/******************************************************/
/*                   randomMapping                    */
/******************************************************/
void randomMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping) {

  //start with initial mapping
  initialMapping(C, P, mapping);

  //randomly exchange positions within mapping[]
  NodeID k = C.number_of_nodes();
  double kdbl = (double) k;
  for(unsigned i = 0; i < k; i++) {
    unsigned newI = (unsigned) (kdbl * rand01());
    if(newI == k) {
      newI--;
    }
    int save = mapping[i];
    mapping[i] = mapping[newI];
    mapping[newI] = save;
  }
}

/*********************** brandfassMapping **************************/
/*                                                                 */
/* Mapping based on what is called "construction method" in        */
/* [Brandfass2012a], and which goes back to [Mueller/Merbach1970a].*/
/* Note that this is only the first part of the mapping method in  */
/* [Brandfass2012a].                                               */
/*                                                                 */
/*******************************************************************/
void brandfassMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping) {

  NodeID k = C.number_of_nodes();
  unsigned int offset = P.number_of_nodes() - k;
  //offset != 0 only for processor graphs with switches

  mapping.resize(k, 0);
  std::vector<int> comFrom(k, 0);
  //total communication involving node i in communication graph
  std::vector<int> procFrom(k, 0);
  //sum of distances to node i in processor graph
  //compute weighted node degrees in communication graph and store them
  //in comFrom
  forall_nodes(C, i) {
    forall_out_edges(C, edge, i) {
      comFrom[i] += C.getEdgeWeight(edge);
    } endfor
	} endfor

	    //compute distances in processor graph
	    //and store them in procFrom
	    for (unsigned int i = 0; i < k; i++) {
	      for (unsigned int j = 0; j < k; j++) {
		procFrom[i] += (P.distanceMatrix)[i + offset][j + offset];
	      }
	    }

  //first assignment: nodeC --> nodeP
  int nodeC = maximumNode(&comFrom);
  int nodeP = minimumNode(&procFrom);
  mapping[nodeC] = nodeP + offset;

  //compute comFrom[i] and procFrom[i],
  //nodeC --> nodeP being the first assignment
  for (unsigned int i = 0; i < k; i++) {
    comFrom[i] = 0;
  }
  forall_out_edges(C, edge, nodeC) {
    comFrom[C.getEdgeTarget(edge)] += C.getEdgeWeight(edge);
  } endfor
      for (unsigned int i = 0; i < k; i++) {
	procFrom[i] = (P.distanceMatrix)[i + offset][nodeP + offset];
      }

  //assigned nodes must not be assigned again
  comFrom[nodeC] = -1;
  procFrom[nodeP] = INT_MAX;

  //Successively assign the not yet assigned node in the communication
  //graph that is most strongly connected to the already assigned
  //nodes to the not yet assigned node in the processor graph that is
  //"most central" w.r.t. the already assigned nodes.

  for (unsigned int i = 1; i < k; i++) {    
    //find new pair (nodeC, nodeP)
    nodeC = maximumNode(&comFrom);
    nodeP = minimumNode(&procFrom);

    //next assignment
    mapping[nodeC] = nodeP + offset;

    if(i < k - 1) {
      //assigned nodes must not be assigned again
      comFrom[nodeC] = -1;
      procFrom[nodeP] = INT_MAX;
      
      //update comFrom
      forall_out_edges(C, edge, nodeC) {
	unsigned int target = C.getEdgeTarget(edge);
	if(comFrom[target] >= 0) {
	  comFrom[target] += C.getEdgeWeight(edge);
	}
      } endfor
	  
	  //update procFrom
	  forall_nodes(C, j) {
	if(procFrom[j] < INT_MAX) {
	  procFrom[j] += (P.distanceMatrix)[j + offset][nodeP + offset];
	}
      } endfor
	  }
  }
}

/********************** brandfassLMapping **************************/
/*                                                                 */
/*             Modification of brandfassMapping:                   */
/* Choose next node in processor graph based on previous choice of */
/* node in communication graph.                                    */
/*                                                                 */
/*******************************************************************/
void brandfassLMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping){
  NodeID k = C.number_of_nodes();
  unsigned int offset = P.number_of_nodes() - k;
  //offset != 0 only for processor graphs with switches
  
  mapping.resize(k, 0);
  std::vector<int> comFrom(k, 0);
  //Total communication involving node i in communication graph
  std::vector<int> procFrom(k, 0);
  //sum of distances to node i in processor graph

  //compute weighted node degrees in communication graph and store them
  //in comFrom
  forall_nodes(C, i) {
    forall_out_edges(C, edge, i) {
      comFrom[i] += C.getEdgeWeight(edge);
    } endfor
	} endfor

	    //compute distances in processor graph
	    //and store them in procFrom
	    for (unsigned int i = 0; i < k; i++) {
	      for (unsigned int j = 0; j < k; j++) {
		procFrom[i] += (P.distanceMatrix)[i + offset][j + offset];
	      }
	    }

  //find new pair (commNode, procNode)
  unsigned int commNode = maximumNode(&comFrom);
  unsigned int procNode = minimumNode(&procFrom);
  
  //reset comFrom[i]
  for (unsigned int i = 0; i < k; i++) {
    comFrom[i] = 0;
  }

  //Successively assign the not yet assigned node in the communication
  //graph that is most strongly connected to the already assigned
  //nodes to the not yet assigned node in the processor graph that is
  //"most central" w.r.t. the already assigned nodes.
  for(unsigned int i = 0; i < k - 1; i++) {    
    //next assignment
    mapping[commNode] = procNode + offset;

    //assigned nodes must not be assigned again
    comFrom[commNode] = -1;
    procFrom[procNode] = INT_MAX;
    
    //update comFrom
    forall_out_edges(C, edge, commNode) {
      unsigned int target = C.getEdgeTarget(edge);
      if(comFrom[target] >= 0) {//i.e., target is unmapped
	comFrom[target] += C.getEdgeWeight(edge);
      }
    } endfor
	
	//find next commNode
	commNode = maximumNode(&comFrom);
    
    //update procFrom
    for (unsigned int j = 1; j < k; j++) {
      if(procFrom[j] < INT_MAX) {
	procFrom[j] = 0;
	forall_out_edges(C, edge, commNode) {
	  unsigned int target = C.getEdgeTarget(edge);
	  if(comFrom[target] < 0) { //i.e., target has been matched already
	    procFrom[j] += (C.getEdgeWeight(edge)) * ((P.distanceMatrix)[j + offset][(mapping[target]) + offset]);
	  }
      	} endfor
	    }
    }
    
    //find next procNode
    procNode = minimumNode(&procFrom);
  }

  //last assignment
  mapping[commNode] = procNode + offset;
}

/************************ hoeflerMapping **************************/
/*                                                                */
/*                    [hoefler2011]-mapping                       */
/*                                                                */
/******************************************************************/
void hoeflerMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping) {

  NodeID k = C.number_of_nodes();
  unsigned int offset = P.number_of_nodes() - k;
  //offset != 0 only for processor graphs in which
  //there are switches
  mapping.resize(k, 0);
  std::vector<int> comFrom(k, 0);
  //total communication involving node i in communication graph
  std::vector<int> procFrom(k, 0);
  //sum of distances to node i in processor graph

  //compute weighted node degrees in communication graph and store them
  //in comFrom
  forall_nodes(C, i) {
    forall_out_edges(C, edge, i) {
      comFrom[i] += C.getEdgeWeight(edge);
    } endfor
	} endfor

	    //choose commNode
	    unsigned int commNode = maximumNode(&comFrom);

  //choose procNode randomly
  unsigned int procNode = (unsigned int) (0.5 + ((double)(k * rand01())) - 1);
  if(procNode == k){
    procNode--;
  };

  //reset comFrom
  for (unsigned int i = 0; i < k; i++) {
    comFrom[i] = 0;
  }

  //Successively assign the not yet assigned node in the communication
  //graph that is most strongly connected to the already assigned
  //nodes to the not yet assigned node in the processor graph that is
  //closest to the already assigned nodes.
  for(unsigned int i = 0; i < k - 1; i++){    
    //next assignment
    mapping[commNode] = procNode + offset;

    //assigned nodes must not be assigned again
    comFrom[commNode] = -1;
    procFrom[procNode] = INT_MAX;
    
    //update comFrom
    forall_out_edges(C, edge, commNode) {
      unsigned int target = C.getEdgeTarget(edge);
      if(comFrom[target] >= 0) {
	if(comFrom[target] < (int) (C.getEdgeWeight(edge))) {
	  comFrom[target] = C.getEdgeWeight(edge);
	}
      }
    } endfor
	
	//find next commNode
	commNode = maximumNode(&comFrom);

    //update procNode: next procNode will be node closest to old procNode
    int minDist = INT_MAX;
    unsigned int newProcNode = procNode;
    forall_nodes(C, j) {
      if(procFrom[j] < INT_MAX) {
	int newDist = (P.distanceMatrix)[j + offset][procNode + offset];
	if(newDist < minDist) {
	  minDist = newDist;
	  newProcNode = j;
	}
      }
    } endfor
	procNode = newProcNode;
  }

  //last assignment
  mapping[commNode] = procNode + offset;

  // std::vector<int> check(k, 0);
  // for(unsigned int i = 0; i < k; i++){
  //   (check[mapping[i]])++;
  // }
  // bool ok = true;
  // for(unsigned int i = 0; i < k; i++){
  //   if(check[mapping[i]] != 1){
  //     ok = false;
  //     printf("mapping[%d] is %d\n", i, mapping[i]);
  //     printf("check[%d] is %d\n", mapping[i], check[mapping[i]]);
  //   }
  // }
  // if(ok == true){
  //   printf("Mapping is valid\n");
  // } else {
  //   printf("INVALID MAPPING\n !!!");
  //   exit(555);
  // }

}

/************************ hoeflerLMapping **************************/
/*                                                                 */
/*             Modification of hoeflerMapping:                     */
/* Choose next node in processor graph based on previous choice of */
/* node in communication graph.                                    */
/*                                                                 */
/*******************************************************************/

void hoeflerLMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping){  
  NodeID k = C.number_of_nodes();
  unsigned int offset = P.number_of_nodes() - k;
  //offset != 0 only for processor graphs in which
  //there are switches
  mapping.resize(k, 0);
  std::vector<int> comFrom(k, 0);
  //total communication involving node i in communication graph
  std::vector<int> procFrom(k, 0);
  //sum of distances to node i in processor graph

  //compute weighted node degrees in communication graph and store them
  //in comFrom
  forall_nodes(C, i) {
    forall_out_edges(C, edge, i) {
      comFrom[i] += C.getEdgeWeight(edge);
    } endfor
	} endfor

	    //choose commNode
	    unsigned int commNode = maximumNode(&comFrom);

  //choose procNode randomly
  unsigned int procNode = (unsigned int) (0.5 + ((double)(k * rand01())) - 1);
  if(procNode == k){
    procNode--;
  };

  //reset comFrom[i]
  for (unsigned int i = 0; i < k; i++) {
    comFrom[i] = 0;
  }
  
  //Successively assign the not yet assigned node in the communication
  //graph that is most strongly connected to the already assigned
  //nodes to the not yet assigned node in the processor graph that is
  //"most central" w.r.t. the already assigned nodes.
  for (unsigned int i = 0; i < k - 1; i++) {    
    //next assignment
    mapping[commNode] = procNode + offset;

    //assigned nodes must not be assigned again
    comFrom[commNode] = -1;
    procFrom[procNode] = INT_MAX;
      
    //update comFrom
    forall_out_edges(C, edge, commNode) {
      unsigned int target = C.getEdgeTarget(edge);
      if(comFrom[target] >= 0) {
	int optionalValue = C.getEdgeWeight(edge);
	if(comFrom[target] < optionalValue) {
	  comFrom[target] = optionalValue;
	}
      }
    } endfor

	//find next commNode
	commNode = maximumNode(&comFrom);
      
    //update procFrom
    for (unsigned int j = 1; j < k; j++) {
      if(procFrom[j] < INT_MAX) {
	procFrom[j] = 0;
	forall_out_edges(C, edge, commNode) {
	  unsigned int target = C.getEdgeTarget(edge);
	  if(comFrom[target] < 0) { //i.e., target has been matched already
	    procFrom[j] += (C.getEdgeWeight(edge)) * ((P.distanceMatrix)[j + offset][(mapping[target]) + offset]);
	  }
	} endfor
	    }
    }
    //find next procNode
    procNode = minimumNode(&procFrom);
  }
  
  //last assignment
  mapping[commNode] = procNode + offset;

  // std::vector<int> check(k, 0);
  // for(unsigned int i = 0; i < k; i++){
  //   (check[mapping[i]])++;
  // }
  // bool ok = true;
  // for(unsigned int i = 0; i < k; i++){
  //   if(check[mapping[i]] != 1){
  //     ok = false;
  //   }
  // }
  // if(ok == true){
  //   printf("Mapping is valid\n");
  // } else {
  //   printf("INVALID MAPPING !!!\n");
  //   exit(555);
}

/************************** almaMapping ****************************/
/*                                                                 */
/*         Amorphous Labelling-based Mapping Algorithm             */
/*                                                                 */
/*******************************************************************/
bool almaMapping(alma_graph_access & g, graph_access & C, processor_graph_access & P,
		 PartitionConfig & map_config, int64_t numberHierarchies,
		 std::vector<int> & mapping){  

  //First, set vlabels[v] to block ID of v and check for consistency
  g.resize_vlabels(g.number_of_nodes());
  int64_t maxPartID = 0;
  forall_nodes(g, v) {
    g.setLabel(v, g.getPartitionIndex(v));
    if(g.getPartitionIndex(v) > maxPartID) {
      maxPartID = g.getPartitionIndex(v);
    }
  } endfor

  vector<int64_t> procVertexLabels;
  bool OK = P.setHammingLabels(P);
  if(OK == false) {
    return(false);
  }

  forall_nodes(g, v) {
    g.setLabel(v, P.getHammingLabel(g.getLabel(v)));
  } endfor

  //At this point, the label of a vertex of g indicates the Hamming label of the vertice's block
  //Labels are now extended to full labels

  g.getVertexLabelsFromPartition(g);
  g.makeVectorLabel2Vertex(g);
  //g.display();
  //g.statisticsOnMissfits();

  //Loop over all hierarchies
  for(int64_t run = 0; run < numberHierarchies; run++) {
    cout << "hey1" << endl;
    Permutation pi(g.getNumDigitsL());
    pi.ini();
    pi.shuffle();
    pi.display();
    Permutation ip(g.getNumDigitsL());
    ip.invert(pi);
        
    cout << "1:::Total communication costs and total cut are " << g.getTotalCommunicationCost(g) <<
      " and " << g.getTotalCut(g) << endl;
    g.permuteVertexLabels(g, pi);
    int64_t digit = g.getNumDigitsL() - 1;//least significant digit
    int64_t level = 0;
    Hierarchy hi(g.getNumDigitsL() - 1);

    cout << "still alive 1xxx" << endl;

    alma_graph_access* newG_ptr = &g;
    while(digit > 0) {// loop over graphs in hierarchy
      cout << "3::newG has " << newG_ptr->number_of_nodes() << " vertices" << endl;
      if(ip.entry(digit) >= g.getNumDigitsForID()) {
	cout << "still alive 1a" << endl;
	newG_ptr->swapVertexLabelsDiversification(*newG_ptr);
	cout << "still alive 1aa" << endl;
      } else {
	//newG->display(*newG_ptr);
	cout << "still alive 1b" << endl;
	newG_ptr->swapVertexLabelsIntensification(*newG_ptr);
	cout << "still alive 1bb" << endl;
      }
      hi.insertGraph(newG_ptr, level);
      vector<NodeID> parents(0);

      cout << "still alive 5" << endl;

      alma_graph_access* h_ptr = new alma_graph_access();
      newG_ptr = newG_ptr->contract(h_ptr, *newG_ptr, parents);

      cout << "still alive 6" << endl;

      cout << "1::newG has " << newG_ptr->number_of_nodes() << " vertices" << endl;

      cout << "still alive 7" << endl;

      hi.insertParents(parents, level);

      cout << "still alive 7.5" << endl;

      cout << "2::newG has " << newG_ptr->number_of_nodes() << " vertices" << endl;

      cout << "still alive 8" << endl;

      digit--;
      level++;
      cout << "2.5::newG has " << newG_ptr->number_of_nodes() << " vertices" << endl;
    }
    cout << "still alive 9" << endl;

    hi.updateLabelsOnFinestLevel();

    cout << "still alive 10" << endl;

    //g = hi.retFinestGraph();
    cout << "11::g has " << g.number_of_nodes() << " vertices" << endl;

    cout << "still alive 11" << endl;

    g.permuteVertexLabels(g, ip);

    cout << "still alive 12" << endl;

    cout << "2:::Total communication costs and total cut are " << g.getTotalCommunicationCost(g) <<
      " and " << g.getTotalCut(g) << endl;
    cout << endl;
  }
  return(true);
}

/*******************************************************************/
/*                     multiSimpleMapping                          */
/*******************************************************************/
void multiSimpleMapping(graph_access & C, processor_graph_access & P, std::vector<int> & finalMapping) {
  
  graph_access coarsedComm;
  C.copy(coarsedComm);

  NodeID k = C.number_of_nodes();
  int level = P.splitDeg.size() - 1;
  bool finished = false;
  int numOfGroups = C.number_of_nodes(); 

  std::vector<std::vector<int> > maps;
  std::vector<std::vector<int> > internalMaps;
  
  while (!finished) {    
  
    if(level == 0) {
      finished = true;
    }

    int degree = P.splitDeg[level];
    numOfGroups /=  degree;

    std::vector<bool> mapped(numOfGroups * degree, false);
    std::vector<int> mapping(numOfGroups * degree, -1);
    std::vector<int> internalPosition(numOfGroups * degree, 0);
    std::queue<int> notFound;


    for (int i = 0; i < numOfGroups; i++) {
      int maxWeight = 0;
      int maxDest = 0;
      int maxTarget = 0;
      
      forall_nodes(coarsedComm, node) {
	if(!mapped[node]) {
	  forall_out_edges(coarsedComm, edge, node) {
	    if(coarsedComm.getEdgeWeight(edge) > maxWeight && (!mapped[coarsedComm.getEdgeTarget(edge)])) {
	      maxWeight = coarsedComm.getEdgeWeight(edge);
	      maxDest = node;
	      maxTarget = coarsedComm.getEdgeTarget(edge);
	    }
	  } endfor
	      }
      } endfor
    
	  if(maxWeight == 0) { 
	    for (int l = 0; l < numOfGroups * degree; l++) {
	      if(!mapped[l]) {
		mapping[l] = i;
		mapped[l] = true;
		for (int k = 1; k < degree; k++) {
		  notFound.push(i);
		}
		break;
	      }
	    }
	  } else {
	    mapping[maxDest] = i;
	    mapped[maxDest] = true;
	    mapping[maxTarget] = i;
	    mapped[maxTarget] = true;
	
	    std::vector<int> connectionStr(numOfGroups * degree, 0);
	
	    forall_out_edges(coarsedComm, edge, maxDest) {
	      connectionStr[coarsedComm.getEdgeTarget(edge)] += coarsedComm.getEdgeWeight(edge);
	    } endfor
		forall_out_edges(coarsedComm, edge, maxTarget) {
	      connectionStr[coarsedComm.getEdgeTarget(edge)] += coarsedComm.getEdgeWeight(edge);
	    } endfor
	
		int size;
	    for (size = 2; size < degree; size++) {
	      int heaviestNode = -1;
	      int weight = -1;
	      for (unsigned int j = 0; j < connectionStr.size(); j++) {
		if(connectionStr[j] > weight && !mapped[j]) {
		  weight = connectionStr[j];
		  heaviestNode = j;
		}
	      }
	      if(weight > 0) {
		mapping[heaviestNode] = i;
		mapped[heaviestNode] = true;
	    
		forall_out_edges(coarsedComm, edge, heaviestNode) {
		  connectionStr[coarsedComm.getEdgeTarget(edge)] += coarsedComm.getEdgeWeight(edge);
		} endfor
		    } else {
		notFound.push(i);
		mapped[i] = true;
	      }
	    }
	  }
    }
    
    std::vector<int> test(numOfGroups, 0);
    for (unsigned int i = 0; i < mapping.size(); i++) {
      if(mapping[i] == -1) {
	mapping[i] = notFound.front();
	notFound.pop();
      }
      test[mapping[i]]++;
    }

    if(P.refine) {
      refinement(coarsedComm, P, mapping, numOfGroups, degree, internalPosition);
    }
    
    for (unsigned int i = 0; i < test.size(); i++) {
      if(test[i] != degree) {
	cout << test[i] << " in Platz " << i << "\n";
	error("Internal Error while mapping.", 150);
      }
    }
    
    maps.push_back(mapping);
    internalMaps.push_back(internalPosition);
    graph_access cc = coarseComm(coarsedComm, mapping, numOfGroups);
   
    cc.copy(coarsedComm);
 
    level--;
  }
  int nextLevel;
  int thisLevel = 1;

  if(maps.size() == 1) {
    nextLevel = k;
  } else {
    nextLevel = k / maps[1].size();
  }

  std::vector<int> oldMap = P.treeSort;
  for (unsigned int i = 0; i < maps.size(); i++) {

    std::vector<int> newMap(k, -1);
    for (unsigned int j = 0; j < k; j++) {
      for (int k = 0; k < nextLevel; k++) {
	if(newMap[maps[i][oldMap[j] / thisLevel] * nextLevel + internalMaps[i][oldMap[j] / thisLevel] * thisLevel + k] == -1) {
	  newMap[maps[i][oldMap[j] / thisLevel] * nextLevel + internalMaps[i][oldMap[j] / thisLevel] * thisLevel +  k] = oldMap[j];
	  break;
	} else {
	  if(k == nextLevel - 1) {
	    cout << internalMaps[i][oldMap[j] / thisLevel] << "\n";
	    cout << j << " ist j." << maps[i][oldMap[j / thisLevel]] * nextLevel << "sollte es hin. \n";
	  }
	}
      }      
    }
    oldMap.swap(newMap);
    thisLevel = nextLevel;
    if(maps.size() <= i + 2) {
      nextLevel = k;
    } else {
      nextLevel = k / maps[i + 2].size();
    }
  }

  //build finalMapping
  finalMapping.swap(oldMap);
  int offset = P.number_of_nodes() - k;
  for (unsigned int i = 0; i < k; i++) {
    finalMapping[i] += offset;
  }
}

/*******************************************************************/
/*                      multiKahipMapping                          */
/*******************************************************************/
void multiKahipMapping(graph_access & C, processor_graph_access & P, const std::string & conf, std::vector<int> & finalMapping) {
  NodeID k = C.number_of_nodes();
  
  graph_access coarsedComm;
  C.copy(coarsedComm);

  int level = P.splitDeg.size() - 1;
  int numOfGroups = k; 
  int nodeNo = k;

  std::vector<std::vector<int> > maps;
  std::vector<std::vector<int> > internalMaps;

  graph_partitioner part;
  PartitionConfig partition_config;
  configuration config;

  if(conf == "fast") {
    config.fast(partition_config);
  } else {
    if(conf == "eco") {
      config.eco(partition_config);
    } else {
      if(conf == "fastsocial") {
	config.fastsocial(partition_config);
      } else {
	if(conf == "ecosocial") {
	  config.ecosocial(partition_config);
	} else {
	  config.fast(partition_config);
	}
      }
    }
  }
 
  partition_config.imbalance = 0;
  partition_config.epsilon = 0;
  partition_config.kaffpa_perfectly_balance = true;
  

  while (level >= 0) {  
    int degree = P.splitDeg[level];

    nodeNo = numOfGroups;
    numOfGroups /= degree;
    
    partition_config.graph_allready_partitioned = false;
    partition_config.k = numOfGroups;  
    coarsedComm.set_partition_count(partition_config.k);
    partition_config.largest_graph_weight = nodeNo;
    partition_config.upper_bound_partition = degree;
    
    part.perform_partitioning(partition_config, coarsedComm);
    
    std::vector<int> mapping;
    std::vector<int> test(numOfGroups, 0);
    std::vector<int> lower;
    std::vector<int> higher;

    forall_nodes(coarsedComm, node) {
      mapping.push_back(coarsedComm.getPartitionIndex(node));   
      test[coarsedComm.getPartitionIndex(node)]++;
      if(test[coarsedComm.getPartitionIndex(node)] > degree) {
	higher.push_back(mapping.size() - 1);
      }
    } endfor

	for (int i = 0; i < numOfGroups; i++) {
	  if(test[i] < degree) {
	    for (int j = test[i];  j < degree; j++) {
	      lower.push_back(i);
	    }
	  }
	}
    
    if(lower.size() != higher.size()) error("Internal error while parsing.",151);
    
    for (unsigned int i = 0; i < lower.size(); i++) {
      mapping[higher[i]] = lower[i];
    }

    std::vector<int> internalPosition(numOfGroups * degree, 0);
    if(P.refine) {
      refinement(coarsedComm, P, mapping, numOfGroups, degree, internalPosition);
    }
    maps.push_back(mapping);
    internalMaps.push_back(internalPosition);
    graph_access cc = coarseComm(coarsedComm, mapping, numOfGroups);
   
    cc.copy(coarsedComm);   
    level--;
  }
  int nextLevel;
  int thisLevel = 1;

  if(maps.size() == 1) {
    nextLevel = k;
  } else {
    nextLevel = k / maps[1].size();
  }
  std::vector<int> oldMap = P.treeSort;
  for (unsigned int i = 0; i < maps.size(); i++) {
    std::vector<int> newMap(k, -1);
    for (unsigned int j = 0; j < k; j++) {
      for (int l = 0; l < nextLevel; l++) {
	if(newMap[maps[i][oldMap[j] / thisLevel] * nextLevel + internalMaps[i][oldMap[j] / thisLevel] * thisLevel +  l] == -1) {
	  newMap[maps[i][oldMap[j] / thisLevel] * nextLevel + internalMaps[i][oldMap[j] / thisLevel] * thisLevel + l] = oldMap[j];
	  break;
	} else {
	  if(l == nextLevel - 1) {
	    cout << "i: " << i << " - " << j << " : j." << maps[i][oldMap[j] / thisLevel] * nextLevel << "sollte es hin. \n";
	  }
	}
      }      
    }
    oldMap.swap(newMap);
    thisLevel = nextLevel;
    if(maps.size() <= i + 2) {
      nextLevel = k;
    } else {
      nextLevel = k / maps[i + 2].size();
    }
  }

  //build final mapping
  finalMapping.swap(oldMap);

  int offset = P.number_of_nodes() - k;
  for (unsigned int i = 0; i < k; i++) {
    finalMapping[i] += offset;
  }
} 
