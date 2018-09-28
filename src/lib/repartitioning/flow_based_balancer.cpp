#include "flow_based_balancer.h"
#include "heap.h"
#include <cstdio>

FlowBased::FlowBased(std::vector<float>** const diff_load_ptr, int max_num_rounds) :
  maxNumPhases(max_num_rounds), diffLoad_ptr(diff_load_ptr) {

}


std::pair<int, float> FlowBased::findBestMove(NodeID vertex, graph_access& G,
				     const std::vector<NodeWeight>& part_weights,
				     const std::vector<NodeWeight>& min_part_weights, 
				     const std::vector<NodeWeight>& max_part_weights,
				     double * const * const pflow) {
  int dest = NO_DEST;
  int v_part = G.getPartitionIndex(vertex);
  float best_prio = NO_PRIO;
  
  // check neighborhood
  forall_out_edges(G,e,vertex) {
    int neigh = G.getEdgeTarget(e);
    int neigh_part = G.getPartitionIndex(neigh);
    float weight = G.getNodeWeight(vertex);
    if (neigh_part != v_part && locked.count(std::make_pair(vertex, neigh_part)) == 0) {
      //debug4("Boundary node %i (weight %f): %i ==> %i\n", vertex, weight, v_part, neigh_part);
      // TODO: test condition on partition sizes!!!
      if (pflow[v_part][neigh_part] > weight * 0.5f ||
	  ((pflow[v_part][neigh_part] > 0.0f) && (part_weights[v_part] > max_part_weights[v_part]) &&
	   part_weights[neigh_part] < min_part_weights[v_part])) {
	float prio = getMovePriorityDiv(vertex, G, v_part, neigh_part);
	if (prio < best_prio) {
	  //debug4("Flow candidate %i: %i ==> %i, prio: %f\n", vertex, v_part, neigh_part, prio);
	  best_prio = prio;
	  dest = neigh_part;
	}
      } else {
	//debug4("%f < %f\n", pflow[v_part][neigh_part], (0.5 * weight));
      }
    }
  } endfor
  
  return std::make_pair(dest, best_prio);
}

void FlowBased::flowBased(PartitionConfig & partition_config, graph_access& G, 
			  const std::vector<NodeWeight>& min_part_weights,
			  const std::vector<NodeWeight>& max_part_weights, 
			  const int numCalls) {

  // initialize part weights and their upper and lower bounds
  int num_parts = partition_config.k;
  int num_nodes = G.number_of_nodes();
  std::vector<NodeWeight> part_weights(num_parts,0);
  forall_nodes(G,n) {
    part_weights[G.getPartitionIndex(n)] += G.getNodeWeight(n);
  } endfor
  int moved_from[num_nodes];
  MemoryFill(moved_from, int, num_nodes);
  
  // determine if balancing necessary
  bool necessary = isBalancingNecessary(min_part_weights, max_part_weights, part_weights);
  if (!necessary) {
    return;// part;
  }
  
  // initialize balancing flow array
  double ** const pflow = MemoryAllocate(double*, num_parts);
  for (int p = 0; p < num_parts; p++)
    pflow[p] = MemoryAllocate(double, num_parts);
  
  // construct priority queue (use heap data structure)
  struct Heap * heap = heap::init(num_nodes);
  
  computeApproximateBalancingFlow(partition_config, G, 
				  min_part_weights, max_part_weights, part_weights, pflow);
  
  int phase = 0;
  while (necessary && phase < maxNumPhases) {

    // populate heap with each eligible node and its priority
    locked.clear();
    for (int v = 0; v < num_nodes; v++) {
      locked.insert(std::make_pair(v, G.getPartitionIndex(v)));
      
      // get priority and destination and insert the data into heap
      std::pair<int, float> dest_prio = findBestMove(v, G, part_weights,
						     min_part_weights, max_part_weights, pflow);

      if (dest_prio.second != NO_PRIO) {
	heap::insert(heap, v, dest_prio.second);
      }
    }
    
    if (heap::size(heap) == 0) {
      //debug3("Heap is empty, quit...\n"); // after displaying moves...\n");
      for (int i = 0; i < num_parts; i++) {
	for (int j = 0; j < num_parts; ++j) {
	  if (pflow[i][j] > 0.5f && numCalls <= 3) {
	    for (int p = 0; p < num_parts; ++p) {
	      MemoryFree(pflow[p], float, num_parts);
	    }
	    MemoryFree(pflow, float*, num_parts);
	    heap::dispose(heap);
	    
	    return flowBased(partition_config, G, min_part_weights, max_part_weights, numCalls + 1);
	  }
	}
      }
      
      for (int p = 0; p < num_parts; ++p) {
	MemoryFree(pflow[p], float, num_parts);
      }
      MemoryFree(pflow, float*, num_parts);
      heap::dispose(heap);
      return;// part;
    }
    
    
    // balancing by node moves in priority order
    while (heap::size(heap) > 0) {
      NodeID vertex = heap::top(heap);
      heap::remove(heap, vertex);
      
      // vertex might change partition
      std::pair<int, float> dest_prio = findBestMove(vertex, G, part_weights, min_part_weights, 
						max_part_weights, pflow);
      
      // check if destination exists
      if (dest_prio.first == NO_DEST) {
	// try to insert element again, hopefully priority has changed since last insert
	if (dest_prio.second != NO_PRIO) {
	  heap::insert(heap, vertex, dest_prio.second);
	}
      } else if (dest_prio.first != moved_from[vertex]) {
	//debug3("balance: vertex %d moved of weight %d from partition %d to partition %d (priority %f)\n", vertex, g.getVertexWeight(vertex), part[vertex], dest_prio.first, heap::priority(heap, vertex));
	//moved_from[vertex] = part[vertex];

	adaptForMigration(G, part_weights, pflow, vertex, dest_prio.first);

	// insert again, seen from new part
	std::pair<int, float> dest_prio = findBestMove(vertex, G, part_weights,
						       min_part_weights, max_part_weights, pflow);
	if (dest_prio.second != NO_PRIO) {
	  heap::insert(heap, vertex, dest_prio.second);
	}
	
	// adapt neighbors
	forall_out_edges(G,e,vertex) {
	  NodeID neigh = G.getEdgeTarget(e);
	  heap::remove(heap, neigh);
	  std::pair<int, float> dest_prio = findBestMove(neigh, G, part_weights, min_part_weights, 
						max_part_weights, pflow);
	  if (dest_prio.second != NO_PRIO) {
	    heap::insert(heap, neigh, dest_prio.second);
	  }
	} endfor
      }
    }

    necessary = false;
    for (int p = 0; p < num_parts; ++p) {
      if (part_weights[p] > max_part_weights[p] || part_weights[p]
	  < min_part_weights[p]) {
	necessary = true;
	break;
      }
    }
    int m = 1;
    for (int p = 0; p < num_parts; ++p) {
      if (part_weights[p] > m) {
	m = part_weights[p];
      }
    }
    //debug3("another balancing step necessary? %i\n", necessary);
    printf("another balancing step necessary? %s\n", necessary ? "true" : "false");
    ++phase;
    printf("phase: %d ... max: %d[%d]\n", phase, m, max_part_weights[0]);
  }
  
  // free memory
  heap::dispose(heap);
  for (int p = 0; p < num_parts; p++) MemoryFree(pflow[p], double, num_parts);
  MemoryFree(pflow, double*, num_parts);
  
  //return part;
}

void FlowBased::adaptForMigration(graph_access& G, std::vector<NodeWeight>& part_weights,
				  double * const * const pflow, NodeID curr_vertex, int dest) {

  //debug4("vertex %i (weight: ", curr_vertex);
  
  // new part found => update flow, part info
  NodeWeight vw = G.getNodeWeight(curr_vertex);

  if (pflow != NULL) {
    //debug4("%i) migrates from %i to %i\n", vw, part[curr_vertex], dest);
    //debug4("old flow: (%i ==> %i): %f\n", part[curr_vertex], dest,
    //	   pflow[part[curr_vertex]][dest]);
    pflow[G.getPartitionIndex(curr_vertex)][dest] -= vw;
    pflow[dest][G.getPartitionIndex(curr_vertex)] += vw;
  }
  part_weights[G.getPartitionIndex(curr_vertex)] -= vw;
  part_weights[dest] += vw;

  
  //debug3("old part[%i]: %i\n", curr_vertex, part[curr_vertex]);
  G.setPartitionIndex(curr_vertex, dest);
  locked.insert(std::make_pair(curr_vertex, dest));
  //debug3("new part[%i]: %i\n", curr_vertex, part[curr_vertex]);
}

void FlowBased::balance(PartitionConfig & partition_config, graph_access& G, 
			const std::vector<NodeWeight>& min_part_weights, 
			const std::vector<NodeWeight>& max_part_weights) {
  /*part = */flowBased(partition_config, G, min_part_weights, max_part_weights,
		   NUM_CALLS_FLOW_BALANCE);
  
  //return part; // balanceSimple(g, min_part_weights, max_part_weights, part);
}

void FlowBased::balanceSimple(PartitionConfig & partition_config, graph_access& G, 
			      const std::vector<NodeWeight>& min_part_weights, 
			      const std::vector<NodeWeight>& max_part_weights) {
  // compute part weights and determine if they need balancing
  std::vector<NodeWeight> part_weights(partition_config.k,0);
  forall_nodes(G,n) {
    part_weights[G.getPartitionIndex(n)] += G.getNodeWeight(n);
  } endfor

  // compute balancing flow
  int P = partition_config.k;
  double ** const pflow = MemoryAllocate(double*, P);
  for (int p = 0; p < P; p++)
    pflow[p] = MemoryClearAllocate(double, P);
  computeApproximateBalancingFlow(partition_config, G, min_part_weights,
				  max_part_weights, part_weights, pflow);


  // construct vector for random order visits
  int V = G.number_of_nodes();
  std::vector<int> visit_order(V);
  for (int v = 0; v < V; ++v) {
    visit_order[v] = v;
  }
  
  int phase = 0;
  while (phase < maxNumPhases && isBalancingNecessary(min_part_weights,
						      max_part_weights, part_weights)) {
    // in each phase: visit vertices in random order
    random_shuffle(visit_order.begin(), visit_order.end());
    for (int v = 0; v < V; ++v) {
      // at each vertex: check (best) possible move(s)
      int curr = visit_order[v];
      std::pair<int, float> dest_prio = findBestMove(curr, G, part_weights,
						     min_part_weights, max_part_weights, pflow);
      
      if (dest_prio.first >= 0) {
	//debug3("RFBalance: Move vertex %i from part %i to %i\n", curr, part[curr], dest_prio.first);
	adaptForMigration(G, part_weights, pflow, curr, dest_prio.first);
      }
    }
    ++phase;
  }

  // free memory
  for (int p = 0; p < P; p++) MemoryFree(pflow[p], double, P);
  MemoryFree(pflow, double*, P);
  
  //return part;
}

bool FlowBased::isBalancingNecessary(const std::vector<NodeWeight>& min_part_weights, 
				     const std::vector<NodeWeight>& max_part_weights, 
				     const std::vector<NodeWeight>& part_weights) const
{
  for (int p = 0; p < static_cast<int> (part_weights.size()); ++p) {
    if (part_weights[p] < min_part_weights[p] || part_weights[p] > max_part_weights[p]) {
      return true;
    }
  }
  return false;
}


void FlowBased::computeApproximateBalancingFlow(PartitionConfig & partition_config, graph_access& G, 
				     const std::vector<NodeWeight>& min_part_weights, 
				     const std::vector<NodeWeight>& max_part_weights, 
				     const std::vector<NodeWeight>& part_weights, 
				     double * const * const pflow) 
{
  int p, v, n, e;
  int P = partition_config.k;
  int V = G.number_of_nodes();
  int ** const cut = MemoryAllocate(int*, P);

  for (p = 0; p < P; p++) {
    cut[p] = MemoryAllocate(int, P);
    MemoryClear(cut[p], int, P);
  }

  for (v = 0; v < V; v++) {
    forall_out_edges(G, e, v) {
      n = G.getEdgeTarget(e);
      if (G.getPartitionIndex(v) != G.getPartitionIndex(n)) {
	cut[G.getPartitionIndex(v)][G.getPartitionIndex(n)]++;
      }
    } endfor
  }

  int pE = 0;
  for (p = 0; p < P; p++)
    for (n = 0; n < P; n++)
      if (cut[p][n] > 0)
	pE++;
  
  double * const pdegree = MemoryAllocate(double, P);
  int * const pof = MemoryAllocate(int, P + 1);
  int * const pto = MemoryAllocate(int, pE);
  double * const pew = MemoryAllocate(double, pE);
  
  // construct partition graph
  e = 0;
  for (p = 0; p < P; p++) {
    pdegree[p] = 0.0f;
    pof[p] = e;
    for (n = 0; n < P; n++)
      if (cut[p][n] > 0) {
	pto[e] = n;
#ifdef BALANCE_UNWEIGHTED
	
	pew[e] = 1.0f;
#else
	
	pew[e] = (double) cut[p][n] * (double) cut[p][n];
#endif
	
	pdegree[p] -= pew[e];
	e++;
      }
  }
  pof[p] = e;
  
  for (p = 0; p < P; p++)
    MemoryFree(cut[p], int, P);
  MemoryFree(cut, int*, P);
  
  double * const lower = MemoryAllocate(double, P);
  double * const upper = MemoryAllocate(double, P);
  
  for (p = 0; p < P; p++) {
    lower[p] = static_cast<double>(ceil(part_weights[p] - min_part_weights[p]));
    upper[p] = static_cast<double>(floor(part_weights[p] - max_part_weights[p]));
  }
  
  double * const f = MemoryAllocate(double, pE);

  HE::solveWeighted(P, pdegree, pof, pto, pew, lower, upper, f, 1e-6);
  
  for (p = 0; p < P; p++)
    MemoryClear(pflow[p], double, P);
  
  for (p = 0; p < P; p++)
    for (e = pof[p]; e < pof[p + 1]; e++)
      pflow[p][pto[e]] = (double) f[e];
  
  MemoryFree(f, double, pE);
  MemoryFree(upper, double, P);
  MemoryFree(lower, double, P);
  MemoryFree(pew, double, pE);
  MemoryFree(pto, int, pE);
  MemoryFree(pof, int, P + 1);
  MemoryFree(pdegree, double, P);
}

/**
 * Solver for weighted balancing flow problem,
 * algorithm based on Hildreth-d'Esopo,
 * cf. Bronstein et al., 3rd ed.,
 * p. 804, Verlag Harri Deutsch, 1997.
 *
 * Minimizes dual problem
 * psi(u) = 1/4 u^T L u + w - \overline{w}
 * iteratively.
 */
void HE::solveWeighted(
		       int V,
		       double * vw,
		       int * of,
		       int * to,
		       double * ew,
		       double * lower,
		       double * upper,
		       double * x,
		       double epsilon)
{
  int i, j, k;

  double * u = MemoryAllocate(double, V);
  double * l = MemoryAllocate(double, V);

  MemoryClear(u, double, V);
  MemoryClear(l, double, V);

  for (k = 0; k < V * V; k++) {
    for (i = 0; i < V; i++) {
      u[i] = -upper[i] - upper[i] - (vw[i] * l[i]);
      for (j = of[i]; j < of[i + 1]; j++)
	u[i] += ew[j] * u[to[j]];
      for (j = of[i]; j < of[i + 1]; j++)
	u[i] -= ew[j] * l[to[j]];
      if (u[i] > 0.0)
	u[i] = 0.0;
      else
	u[i] /= -vw[i];
    }
    for (i = 0; i < V; i++) {
      l[i] = lower[i] + lower[i] + (vw[i] * u[i]);
      for (j = of[i]; j < of[i + 1]; j++)
	l[i] -= ew[j] * u[to[j]];
      for (j = of[i]; j < of[i + 1]; j++)
	l[i] += ew[j] * l[to[j]];
      if (l[i] > 0.0)
	l[i] = 0.0;
      else
	l[i] /= -vw[i];
    }
  }

  for (i = 0; i < V; i++) {
    for (j = of[i]; j < of[i + 1]; j++)
      x[j] = -0.5 * ew[j] * (u[i] - l[i] + (l[to[j]] - u[to[j]]));
  }

  MemoryFree(l, double, V);
  MemoryFree(u, double, V);
}
