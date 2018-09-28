#ifndef FLOW_BASED_BALANCER_H_
#define FLOW_BASED_BALANCER_H_

#include <algorithm>
#include <set>
#include <cfloat>
#include <vector>
#include <cstring>

#include "data_structure/graph_access.h"
#include "partition_config.h"

#define MemoryAllocate(T, n) (T *)malloc((n) * sizeof(T))
#define MemoryRelocate(p, T, n) (T *)realloc((void *)p, ((n) * sizeof(T)))
#define MemoryFree(p, T, n) free((void *)p)
#define MemoryCopy(d, s, T, n) memcpy((void *)d, (void *)s, (n) * sizeof(T))
#define MemoryClear(p, T, n) memset((void *)p, 0, (n) * sizeof(T))
#define MemoryFill(p, T, n) memset((void *)p, -1, (n) * sizeof(T))
#define MemoryClearAllocate(T, n) (T *)calloc((n), sizeof(T))
#define SafeMemoryAllocate(p, T, n) { assert((p = MemoryAllocate(T, n)) != NULL); }
#define SafeMemoryFree(p, T, n) {if(p) {free((void *) p); (p)=NULL;}}
#define SafeMemoryClearAllocate(p, T, n) { assert((p = MemoryClearAllocate(T, n)) != NULL); }
#define SafeMemoryDelete(p) { if(p) { delete(p); (p)=NULL; }}

#define NUM_CALLS_FLOW_BALANCE 1
#define NUM_RANDOM_PHASES 50

const int NO_DEST = -1;
const float NO_PRIO = FLT_MAX;

class FlowBased {
private:
  int maxNumPhases; ///> Maximum number of iterations for the balancing process.
  std::vector<float>** diffLoad_ptr; ///> Array of pointers to load vectors.
  std::set<std::pair<NodeID, int> > locked;
  
  void balanceGreedy(const graph_access& G, const std::vector<NodeWeight>& min_part_weights,
		     const std::vector<NodeWeight>& max_part_weights);

protected:
  
  inline float getMovePriority(NodeID vertex, const graph_access& G, int src, int dest) {
    return (*diffLoad_ptr[dest])[vertex] - (*diffLoad_ptr[src])[vertex];
  }

  inline float getMovePriorityDiv(NodeID vertex, const graph_access& G, int src, int dest) {
    return (*diffLoad_ptr[dest])[vertex] / (*diffLoad_ptr[src])[vertex];
  }
	
  inline float getMovePriorityAbs(NodeID vertex, const graph_access& G, int src, int dest) {
    return (*diffLoad_ptr[dest])[vertex];
  }

  std::pair<int, float> findBestMove(NodeID vertex, graph_access& G,
				     const std::vector<NodeWeight>& part_weights,
				     const std::vector<NodeWeight>& min_part_weights, 
				     const std::vector<NodeWeight>& max_part_weights,
				     double * const * const pflow);

  void balanceSimple(PartitionConfig & partition_config, graph_access& G, 
		     const std::vector<NodeWeight>& min_part_weights, 
		     const std::vector<NodeWeight>& max_part_weights);

  void flowBased(PartitionConfig & partition_config, graph_access& G, 
		 const std::vector<NodeWeight>& min_part_weights,
		 const std::vector<NodeWeight>& max_part_weights, 
		 const int numCalls = NUM_CALLS_FLOW_BALANCE);

  void adaptForMigration(graph_access& G, std::vector<NodeWeight>& part_weights,
			 double * const * const pflow, NodeID curr_vertex, int dest);

public:

  FlowBased(std::vector<float>** const diff_load_ptr, int max_num_rounds = NUM_RANDOM_PHASES);

  void balance(PartitionConfig & partition_config, graph_access& G, 
	       const std::vector<NodeWeight>& min_part_weights, 
	       const std::vector<NodeWeight>& max_part_weights);

  bool isBalancingNecessary(const std::vector<NodeWeight>& min_part_weights, 
				       const std::vector<NodeWeight>& max_part_weights, 
				       const std::vector<NodeWeight>& part_weights) const;

  void computeApproximateBalancingFlow(PartitionConfig & partition_config, graph_access& G, 
				       const std::vector<NodeWeight>& min_part_weights, 
				       const std::vector<NodeWeight>& max_part_weights, 
				       const std::vector<NodeWeight>& part_weights, 
				       double * const * const pflow);

};

class HE {

public:
  static void solveWeighted(int V,
			    double * vw,
			    int * of,
			    int * to,
			    double * ew,
			    double * lower,
			    double * upper,
			    double * x,
			    double epsilon);

};


#endif /* FLOW_BASED_BALANCER_H_ */
