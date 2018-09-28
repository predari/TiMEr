//code based upon dibap-lite's trunccons.cpp and scale.h

//TODO: FlowBased
//      Smooth

#include "diffusive_refinement.h"
#include "partition/coarsening/clustering/node_ordering.h"
#include "tools/random_functions.h"
#include "quality_metrics.h"
#include "flow_based_balancer.h"
#include <vector>
#include <set>
#include <queue>
#include <algorithm>
#include <cmath>

diffusive_refinement::diffusive_refinement() {
  /*numConsols = 9;
    numDiffIters = 14;	*/
  useEdgeWeights = true;
  useEdgeWeightsInSmooth = true;
  /*numScaleIters = 16;
  useAllNeighbors = true;
  sizeOfNeighborhood = 8;
  limited = true;
  scaleBalancing = true;
  flowbasedBalancing = false;
  greedyBalancing = true;*/
}

diffusive_refinement::~diffusive_refinement() {
                
}

EdgeWeight diffusive_refinement::perform_refinement(PartitionConfig & repartition_config, 
						    graph_access & G, 
						    complete_boundary & boundary) {
  trunccons(repartition_config, G);
  return 0; //TODO return improvement?
}

void diffusive_refinement::trunccons(PartitionConfig & repartition_config, graph_access & G) {
  NodeID num_nodes = G.number_of_nodes();
  int maxDegree = G.getMaxDegree();
  int num_parts = repartition_config.k;
  std::vector<int> partSizes(num_parts);
  std::vector<float>** w; // pointer to the load-vectors for each partition
  std::vector<float>** temp; // pointer to do calculations on the load.
  //active nodes:
  std::vector<bool>** active;
  std::vector<bool>** newActive;

  w = (std::vector<float>**)malloc((num_parts) * sizeof(std::vector<float>*));
  temp = (std::vector<float>**)malloc((num_parts) * sizeof(std::vector<float>*));
  active = (std::vector<bool>**)malloc((num_parts) * sizeof(std::vector<bool>*));
  newActive = (std::vector<bool>**)malloc((num_parts) * sizeof(std::vector<bool>*));

  for (int c = 0; c < num_parts; c++) {
    w[c] = new std::vector<float> (num_nodes);
    temp[c] = new std::vector<float> (num_nodes);
    active[c] = new std::vector<bool> (num_nodes);
    newActive[c] = new std::vector<bool> (num_nodes);
  }

  std::vector<std::set<PartitionID> > neighborhood(num_parts);
  if (repartition_config.trunccons_limited) {
    std::vector<std::vector<std::pair<int,PartitionID> > > count(num_parts);
    for (int i = 0; i < num_parts; ++i) {
      for (int j = 0; j < num_parts; ++j) {
	count[i].push_back(std::pair<int,PartitionID>(0,j));
      }
    }
    forall_nodes(G, node) {
      count[G.getPartitionIndex(node)][G.getPartitionIndex(node)].first++;
      forall_out_edges(G, e, node) {
	count[G.getPartitionIndex(node)][G.getPartitionIndex(G.getEdgeTarget(e))].first++;
      } endfor
    } endfor
    for (int i = 0; i < num_parts; ++i) {
      std::sort(count[i].begin(), count[i].end());
      int m = repartition_config.trunccons_limited_useAllNeighbors ? num_parts 
	: repartition_config.trunccons_limited_sizeOfNeighborhood;
      for (int j = 0; j < m; ++j) {
	int x = count[i].size() - 1 - j;
	if (count[i][x].first == 0) break;
	neighborhood[i].insert(count[i][x].second);
      }
    }
  }
  float alpha = 1.0 / (float) (maxDegree + 1);

  for (int lambda = 0; lambda < repartition_config.trunccons_numConsols; lambda++) {
    bool changed = false;
    PartitionID p;
    float m;

    //calculate partsizes
    for (int i = 0; i < num_parts; ++i) partSizes[i] = 0;
    forall_nodes(G, node) {
      partSizes[G.getPartitionIndex(node)] += 1;
    } endfor

    for (int c = 0; c < num_parts; c++) {
      initLoad(num_nodes, partSizes[c], w[c], c, active[c], G);
    }

    //for all existent partitions
    for (int c = 0; c < num_parts; c++) {
      loadExchange(repartition_config, useEdgeWeights, w[c], temp[c], num_nodes, alpha, active[c], newActive[c], G, neighborhood[c]);
    }

    //the load exchange is done.
    //Find out from which partition each node got the highest load from:
    forall_nodes (G, v) {
      PartitionID oldp = G.getPartitionIndex(v);
      p = 0;
      std::vector<float>& deRef1 = *w[0];
      m = deRef1[v];
      if (repartition_config.trunccons_limited) {
	for (std::set<PartitionID>::iterator it = neighborhood[oldp].begin(); 
	     it  != neighborhood[oldp].end(); ++it) {
	  int c = *it;
	  std::vector<float>& deRef = *w[c];
	  if (deRef[v] > m) {
	    m = deRef[v];
	    p = c;
	  }
	}
      } 
      else {
	for (int c = 1; c < num_parts; c++) {
	  std::vector<float>& deRef = *w[c];
	  if (deRef[v] > m) {
	    m = deRef[v];
	    p = c;
	  }
	}
      }
      //update the partition information for node v.
      if (G.getPartitionIndex(v) != p) {
	G.setPartitionIndex(v, p);
	changed = true;
	forall_out_edges(G, e, v) {
	  //neighborhood[G.getPartitionIndex(v)].insert(G.getPartitionIndex(G.getEdgeTarget(e)));
	} endfor
	    // TODO: erase old neighbors? Caution: don't erase too much
      }
    } endfor

    //end trunccons -- no node has been changed.
    if (!changed) {
      break;
    }
    if (repartition_config.scale_balancing && lambda % BAL_ROUNDS == 0)
      scalebalancing(w, repartition_config, G);

  } //for lambda
  if (repartition_config.scale_balancing) scalebalancing(w, repartition_config, G);
  if (repartition_config.flow_based_balancing) {
    FlowBased flow_based(w, 50);
    std::vector<NodeWeight> max_part_weights(num_parts);
    std::vector<NodeWeight> min_part_weights(num_parts);
    for (int i = 0; i < num_parts; ++i) {
      max_part_weights[i] = repartition_config.upper_bound_partition;
      min_part_weights[i] = 1; //TODO: ?
    }
    flow_based.balance(repartition_config, G, min_part_weights, max_part_weights);
  }
  if (repartition_config.greedy_balancing) greedy_balancing(repartition_config, G, w);
  /*
  // smooth boundaries
  Smooth<float> smoother(w, useEdgeWeightsInSmooth);
  part = smoother.smooth(g, part);
  if (levelNumber == 0) {
    part = flow_based.balance(g, minPartWeights, maxPartWeights, part);
  }
  */

  for (unsigned int c = 0; c < repartition_config.k; c++) {
    delete(w[c]);
    delete(temp[c]);
    delete(active[c]);
    delete(newActive[c]); 
  }
  free((void *)temp);
  free((void *)w);
  free((void *)active);
  free((void *)newActive);
}

void diffusive_refinement::loadExchange(PartitionConfig& repartition_config,
					bool useEdgeWeights,
					std::vector<float>* load,
					std::vector<float>* newLoad,
					const int num_nodes,
					const double alpha,
					std::vector<bool>* active,
					std::vector<bool>* newActive,
					graph_access & G,
					std::set<PartitionID>& neighborhood)
{
  //load is the initial load-vector.
  //newLoad is used as a tempoaray variable
  std::vector<float>* helpSwap;
  float flow;
  EdgeWeight degree;

  //	EdgeWeight* ew = g->getEdgeWeights();
  EdgeWeight* ew = (EdgeWeight*)malloc(sizeof(EdgeWeight) *G.number_of_edges());
  forall_edges(G, e) {
    ew[e] = G.getEdgeWeight(e);
  } endfor

  std::vector<bool>* helpActiveSwap;

  for (int t = 0; t < repartition_config.trunccons_numDiffIters; t++) {
    std::vector<bool>& newlyActive = *newActive;
    std::vector<bool>& active_ref = *active;

    //all nodes that are active now are active afterwards...
    for (int v = 0; v < num_nodes; v++) {
      newlyActive[v] = active_ref[v];
    }

    std::vector<float>& load_ref = *load;
    std::vector<float>& newLoad_ref = *newLoad;

    if (!useEdgeWeights || ew == NULL) {
      //load exchange without edge weights
      //for (int v = 0; v < num_nodes; v++) {
      forall_nodes(G,v) {
	if (active_ref[v] == true && 
	    (!repartition_config.trunccons_limited 
	     || neighborhood.find(G.getPartitionIndex(v)) != neighborhood.end())) {
	  flow = G.getNodeDegree(v) * load_ref[v];
	  forall_out_edges(G,e,v) {
	    flow -= load_ref[G.getEdgeTarget(e)];
	    newlyActive[G.getEdgeTarget(e)] = true;
	  } endfor
	  newLoad_ref[v] = load_ref[v] - alpha * flow;
	}
	else {
	  newLoad_ref[v] = load_ref[v];
	}
      } endfor

    }
    else {

      //load exchange with edge weights
      for (int v = 0; v < num_nodes; v++) {
	if (active_ref[v] == true && 
	    (!repartition_config.trunccons_limited 
	     || neighborhood.find(G.getPartitionIndex(v)) != neighborhood.end())) {
	  degree = 0;
	  flow = 0.0;
	  //calculating the load arriving at node v:
	  forall_out_edges(G,e,v) {
	    degree += G.getEdgeWeight(e);
	    flow -= load_ref[G.getEdgeTarget(e)] * ((float) G.getEdgeWeight(e));

	    newlyActive[G.getEdgeTarget(e)] = true;
	  } endfor
	      
	  flow += (float) degree * load_ref[v];
	  flow *= alpha;
	  newLoad_ref[v] = load_ref[v] - flow;
	}
      
	else {
	  newLoad_ref[v] = load_ref[v];
	}
      }
    }

    //now temp is the pointer to the new calculated load.
    helpSwap = load;
    load = newLoad;
    newLoad = helpSwap;

    helpActiveSwap = active;
    active = newActive;
    newActive = helpActiveSwap;
  }

  free((void *)ew);
}

void diffusive_refinement::initLoad(
	      int num_nodes,
	      int partSize,
	      std::vector<float>* load,
	      PartitionID partitionNumber,
	      std::vector<bool>* active,
	      graph_access & G)
{
  //set the initial loads
  float num_nodes_float = (float) num_nodes;
  std::vector<float>& deRef = *load;

  std::vector<bool>& active_ref = *active;

  //init all nodes to be inactive
  for (int v = 0; v < num_nodes; v++) {
    active_ref[v] = false;
  }

  for (int v = 0; v < num_nodes; v++) {
    if (G.getPartitionIndex(v) == partitionNumber) {
      deRef[v] = num_nodes_float / (float) partSize;
      //iterate over all neighbors of v:
      forall_out_edges(G,e,v) {
	if (G.getPartitionIndex(G.getEdgeTarget(e)) != partitionNumber) {
	  //to[u] is a neighbor in another partition - thus v must be active!
	  active_ref[v] = true;
	  active_ref[G.getEdgeTarget(e)] = true;
	}
      } endfor
    }
    else {
      deRef[v] = 0.0;
    }
  }

}

bool diffusive_refinement::scalebalancing(std::vector<float>** diffLoad_ptr,
					  PartitionConfig & repartition_config, 
					  graph_access & G) {
  int num_parts = repartition_config.k;
  std::vector<NodeWeight> part_weights(num_parts);
  std::vector<float> scales;
  std::vector<float> old_scales;
  //  Partition tmp_part(part);
  std::vector<NodeID> tmp_part(G.number_of_nodes());
  forall_nodes(G, n) {
    tmp_part[n] = G.getPartitionIndex(n);
  } endfor

  std::vector<NodeWeight> max_part_weights(num_parts);
  for (int i = 0; i < num_parts; ++i) 
    max_part_weights[i] = repartition_config.upper_bound_partition;

  // initialize all scales with 1
  scales.assign(num_parts, 1.0);

  // "learn" scales for balanced parts
  for (int i = 0; i < repartition_config.scale_balancing_iters; ++i) {
    old_scales = scales;
    //part_weights = tmp_part.getPartWeights(g, part_weights);
    part_weights.assign(num_parts, 0);
    // determine part weights
    for (unsigned int j = 0; j < tmp_part.size(); ++j) {
      part_weights[tmp_part[j]] += G.getNodeWeight(j);
    }

    /*if (!isBalancingNecessary(min_part_weights, max_part_weights, part_weights)) {
      debug3("Scale balancing finished early in iteration %i.\n", i);
      break;
    }*/

    // compute a new scale for each part
    scales = updateScales(part_weights, max_part_weights, scales);

    // update partition
    tmp_part = scaleBasedPartition(diffLoad_ptr, scales, tmp_part);
    //		scales = normalizeFactors(scales);

    bool hasEmptyParts = false;
    std::vector<bool> emptyParts(num_parts, true);
    for (unsigned int j = 0; j < tmp_part.size(); ++j) {
      emptyParts[tmp_part[j]] = false;
    }
    for (unsigned int j = 0; j < emptyParts.size(); ++j) {
      if (emptyParts[j]) { hasEmptyParts = true; break; }
    }

    if (hasEmptyParts) {
      return false;
    }
  }

  quality_metrics qm;
  if (imbalance(G, tmp_part) < qm.balance(G)) {
    adaptLoads(diffLoad_ptr, scales);
    forall_nodes(G, n) {
      //if (G.getPartitionIndex(n) != tmp_part[n]) std::cout << "scale based change\n";
      G.setPartitionIndex(n, tmp_part[n]);
    } endfor
  }
  
  return true;
}

float diffusive_refinement::imbalance(graph_access& G, std::vector<NodeWeight> part) {
  std::vector<PartitionID> part_weights(part.size(), 0);

  double overallWeight = 0;

  forall_nodes(G, n) {
    PartitionID curPartition = part[n];
    part_weights[curPartition] += G.getNodeWeight(n);
    overallWeight += G.getNodeWeight(n);
  } endfor

  double balance_part_weight = ceil(overallWeight / (double)G.get_partition_count());
  double cur_max             = -1;

  forall_blocks(G, p) {
    double cur = part_weights[p];
    if (cur > cur_max) {
      cur_max = cur;
    }
  } endfor

  double percentage = cur_max/balance_part_weight;
  return percentage;
}

std::vector<float> diffusive_refinement::updateScales(
                                         std::vector<NodeWeight>& part_weights, 
					 std::vector<NodeWeight>& max_part_weights, 
					 std::vector<float>& scales) const
{
  float q = 0.0f;
  float q2 = 0.0f;
  const float one_minus_damp = 1.0 - SCALE_BALANCE_DAMP;

  for (int p = 0; p < (int) scales.size(); p++) {
    //debug4("max[%i]: %i, curr[%i]: %i\n", p, max_part_weights[p], p, part_weights[p]);

    q = static_cast<float> (max_part_weights[p]) / static_cast<float> (part_weights[p]);
    q2 = q*q;

    if (q2 > 1.667f)
      q2 = 1.667f;

    if (q2 < 0.6f)
      q2 = 0.6f;

    //    scales[p] *= SCALE_BALANCE_DAMP + one_minus_damp * q2;
    scales[p] = scales[p] * SCALE_BALANCE_DAMP + one_minus_damp * q2;

    if (scales[p] > SCALE_BALANCE_LIMIT)
      scales[p] = SCALE_BALANCE_LIMIT;

    if (scales[p] < 1.0f / SCALE_BALANCE_LIMIT)
      scales[p] = 1.0f / SCALE_BALANCE_LIMIT;
  }

  return scales;
}

std::vector<NodeID> diffusive_refinement::scaleBasedPartition(std::vector<float>** diffLoad_ptr,
							      std::vector<float>& scales,
							      std::vector<NodeID> partition) {
  double maxi = 0.0;
  int argmax = 0;

  // TODO: rethink loop order, not cache-efficient right now!
  for (int v = 0; v < (int) partition.size(); ++v) {
    argmax = 0;
    maxi = scales[0] * (*diffLoad_ptr[0])[v];
    for (int p = 1; p < (int) scales.size(); ++p) {
      if (scales[p] * (*diffLoad_ptr[p])[v] > maxi) {
	argmax = p;
	maxi = scales[p] * (*diffLoad_ptr[p])[v];
      }
    } 
    partition[v] = argmax;
  }
  return partition;
}

void diffusive_refinement::adaptLoads(std::vector<float>** diffLoad_ptr, 
				      std::vector<float>& scales) {
  for (int p = 0; p < (int) scales.size(); ++p) {
    float scale_val = scales[p];
    for (int v = 0; v < (int) diffLoad_ptr[p]->size(); ++v) {
      (*diffLoad_ptr[p])[v] *= scale_val;
    }
  }
}

void diffusive_refinement::greedy_balancing(PartitionConfig& repartition_config,
					    graph_access& G,
					    std::vector<float>** diffLoad_ptr) {
  std::vector<NodeWeight> part_weights(repartition_config.k, 0);
  forall_nodes(G,n) {
    part_weights[G.getPartitionIndex(n)] += G.getNodeWeight(n);
  } endfor

  NodeWeight maxpartweight = repartition_config.upper_bound_partition;
  const PartitionID NOPART = repartition_config.k;
  for (int i = 0; i < 1; ++i) {
    std::priority_queue<std::pair<float,std::pair<NodeID,PartitionID> > > queue;
    forall_nodes(G,n) {
      PartitionID part = G.getPartitionIndex(n);
      if (part_weights[part] > maxpartweight) {
	float maxp = -1;
	PartitionID best = NOPART;
	forall_out_edges(G,e,n) {
	  PartitionID c = G.getPartitionIndex(G.getEdgeTarget(e));
	  if (c == part) continue;
	  float p = (part_weights[part]-part_weights[c]);
	  p= p * (*diffLoad_ptr[c])[n];
	  if (maxp < p) {
	    maxp = p;
	    best = c;
	  }
	} endfor
	if (best != NOPART) {
	  queue.push(std::make_pair(maxp, std::make_pair(n,best)));
	}
      }
    } endfor

    while (!queue.empty()) {
      NodeID n = queue.top().second.first;
      PartitionID best = queue.top().second.second;
      PartitionID part = G.getPartitionIndex(n);
      if (part_weights[part] > maxpartweight) {
	part_weights[part] -= G.getNodeWeight(n);
	part_weights[best] += G.getNodeWeight(n);
	G.setPartitionIndex(n, best);
      }
      queue.pop();
    }
  }
}
