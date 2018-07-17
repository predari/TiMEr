#ifndef DIFFUSIVE_REFINEMENT_H
#define DIFFUSIVE_REFINEMENT_H

#include "uncoarsening/refinement/refinement.h"
#include "data_structure/graph_access.h"
#include "../karma_config.h"
#include <set>
#include <vector>

const int BAL_ROUNDS = 4;
const float SCALE_BALANCE_DAMP = 0.9;
const float SCALE_BALANCE_LIMIT = 4.0;

class diffusive_refinement : public refinement {
private:
  bool useEdgeWeights;	        // use edgeweights in trunccons
  bool useEdgeWeightsInSmooth;	// use edgeweights in smoother.

  //PartitionConfig doCast(PartitionConfig & partition_config);

  void initLoad(int num_nodes, int partSize, std::vector<float>* load, PartitionID partitionNumber,
		std::vector<bool>* active, graph_access & G);

  bool scalebalancing(std::vector<float>** diffLoad_ptr,
		      PartitionConfig & partition_config, 
		      graph_access & G);
  std::vector<float> updateScales(std::vector<NodeWeight>& part_weights, 
				  std::vector<NodeWeight>& max_part_weights, 
				  std::vector<float>& scales) const;
  std::vector<NodeID> scaleBasedPartition(std::vector<float>** diffLoad_ptr,
					  std::vector<float>& scales,
					  std::vector<NodeID> partition);
  void adaptLoads(std::vector<float>** diffLoad_ptr, 
		  std::vector<float>& scales);
  float imbalance(graph_access & G,std::vector<NodeWeight> part);
  
  void trunccons(PartitionConfig & partition_config, 
			graph_access & G);
  void loadExchange(PartitionConfig& config,
		      bool useEdgeWeights,
		      std::vector<float>* load,
		      std::vector<float>* newLoad,
		      const int num_nodes,
		      const double alpha,
		      std::vector<bool>* active,
		      std::vector<bool>* newActive,
		      graph_access & G,
		      std::set<PartitionID>& neighborhood);

  void greedy_balancing(PartitionConfig& partition_config,
				   graph_access& G,
				   std::vector<float>** diffLoad_ptr);

public:
  diffusive_refinement();
  virtual ~diffusive_refinement();

  virtual EdgeWeight perform_refinement(PartitionConfig & config, 
					graph_access& G, 
					complete_boundary & boundary); 
};

#endif //include-guard DIFFUSIVE_REFINEMENT_H
