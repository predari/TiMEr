//
// Author: Christian Schulz <christian.schulz@kit.edu>
// Manipulator: Roland Glantz (RG)// 

#include <algorithm>
#include <cmath>

#include "extended_quality_metrics.h"

extended_quality_metrics::extended_quality_metrics() {
}

extended_quality_metrics::~extended_quality_metrics () {
}

EdgeWeight extended_quality_metrics::edge_cut(graph_access & G) {
        EdgeWeight edgeCut = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = G.getPartitionIndex(n);
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);

                        if (partitionIDSource != partitionIDTarget) {
                                edgeCut += G.getEdgeWeight(e);
                        }
                } endfor 
        } endfor
        return edgeCut/2;
}

EdgeWeight extended_quality_metrics::edge_cut(graph_access & G, int * partition_map) {
        EdgeWeight edgeCut = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = partition_map[n];
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = partition_map[targetNode];

                        if (partitionIDSource != partitionIDTarget) {
                                edgeCut += G.getEdgeWeight(e);
                        }
                } endfor 
        } endfor
        return edgeCut/2;
}

EdgeWeight extended_quality_metrics::edge_cut(graph_access & G, PartitionID lhs, PartitionID rhs) {
        EdgeWeight edgeCut = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = G.getPartitionIndex(n);
                if(partitionIDSource != lhs) continue;
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);

                        if(partitionIDTarget == rhs) {
                                edgeCut += G.getEdgeWeight(e);
                        }
                } endfor 
        } endfor
        return edgeCut;
}

EdgeWeight extended_quality_metrics::edgeCuts(graph_access & G, std::vector<EdgeWeight> & block_internal_edges, std::vector<EdgeWeight> & block_cut) {
  EdgeWeight edgeCut = 0;
  forall_nodes(G, n){ 
    PartitionID partitionIDSource=G.getPartitionIndex(n);
    forall_out_edges(G, e, n){
      NodeID targetNode = G.getEdgeTarget(e);
      PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);
      if (partitionIDSource != partitionIDTarget){
	edgeCut++;
	(block_cut[partitionIDSource])++;
	// edgeCut+= G.getEdgeWeight(e);
	// block_cut[partitionIDSource]+=G.getEdgeWeight(e);
      }
      else{
	if(n < targetNode){
	  block_internal_edges[partitionIDSource]++;
	}
      }
    } endfor 
  } endfor
  return edgeCut/2;
}

void extended_quality_metrics::invade(graph_access & G, NodeID n, PartitionID pID, std::vector<bool> & visited){
  visited[n] = true;
  forall_out_edges(G, e, n){
    NodeID t = G.getEdgeTarget(e);
    if((visited[t] == false) && (G.getPartitionIndex(t) == pID)){
      invade(G, t, pID, visited);
    }
  } endfor
}


NodeID extended_quality_metrics::blockComps(graph_access & G) {
  NodeID nrComps = 0;
  std::vector<bool> visited(G.number_of_nodes(), false);
  forall_nodes(G, n){ 
    if(visited[n] == false){
      PartitionID pID=G.getPartitionIndex(n);
      invade(G, n, pID, visited);
      nrComps ++;
    }
  } endfor
  return nrComps;
}

EdgeWeight extended_quality_metrics::max_communication_volume(graph_access & G, int * partition_map) {
    std::vector<EdgeWeight> block_volume(G.get_partition_count(),0);
    forall_nodes(G, node) {
        PartitionID block = partition_map[node];
        std::vector<bool> block_incident(G.get_partition_count(), false);
        block_incident[block] = true;

        int num_incident_blocks = 0;

        forall_out_edges(G, e, node) {
            NodeID target = G.getEdgeTarget(e);
            PartitionID target_block = partition_map[target];
            if(!block_incident[target_block]) {
                block_incident[target_block] = true;
                num_incident_blocks+=G.getNodeWeight(node);
            }
        } endfor
        block_volume[block] += num_incident_blocks;
    } endfor

    EdgeWeight max_comm_volume = *(std::max_element(block_volume.begin(), block_volume.end()));
    return max_comm_volume;
}

EdgeWeight extended_quality_metrics::max_communication_volume(graph_access & G){
  std::vector<EdgeWeight> block_volume(G.get_partition_count(),0);
  forall_nodes(G, node) {
    PartitionID block = G.getPartitionIndex(node);
    std::vector<bool> block_incident(G.get_partition_count(), false);
    block_incident[block] = true;
    int num_incident_blocks = 0;
    
    forall_out_edges(G, e, node) {
      NodeID target = G.getEdgeTarget(e);
      PartitionID target_block = G.getPartitionIndex(target);
      if(!block_incident[target_block]) {
	block_incident[target_block] = true;
	num_incident_blocks+=G.getNodeWeight(node);
      }
    } endfor
    block_volume[block] += num_incident_blocks;
  } endfor
      
  EdgeWeight max_comm_volume = *(std::max_element(block_volume.begin(), block_volume.end()));
  return max_comm_volume;
}

EdgeWeight extended_quality_metrics::communicationVolumes(graph_access & G, std::vector<EdgeWeight> & block_volume) {
  forall_nodes(G, node) {
    PartitionID block = G.getPartitionIndex(node);
    std::vector<bool> block_incident(G.get_partition_count(), false);
    block_incident[block] = true;
    int num_incident_blocks = 0;
    
    forall_out_edges(G, e, node) {
      NodeID target = G.getEdgeTarget(e);
      PartitionID target_block = G.getPartitionIndex(target);
      if(!block_incident[target_block]) {
	block_incident[target_block] = true;
	num_incident_blocks++;
      }
    } endfor
    block_volume[block] += num_incident_blocks;
  } endfor
      
  EdgeWeight max_comm_volume = *(std::max_element(block_volume.begin(), block_volume.end()));
  return max_comm_volume;
}


int extended_quality_metrics::boundary_nodes(graph_access& G){
  int no_of_boundary_nodes = 0;
  forall_nodes(G, n) { 
    PartitionID partitionIDSource = G.getPartitionIndex(n);
    
    forall_out_edges(G, e, n) {
      NodeID targetNode = G.getEdgeTarget(e);
      PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);
      
      if (partitionIDSource != partitionIDTarget) {
	no_of_boundary_nodes++;
	break; 
      }
    } endfor 
 } endfor
 return no_of_boundary_nodes;
}

int extended_quality_metrics::boundaryNodes(graph_access & G, std::vector<NodeID> & block_nodes, std::vector<NodeID> & block_bnd){
  int no_of_boundary_nodes = 0;
  forall_nodes(G, n) { 
    PartitionID partitionIDSource = G.getPartitionIndex(n);
    block_nodes[partitionIDSource]++;   
    forall_out_edges(G, e, n) {
      NodeID targetNode = G.getEdgeTarget(e);
      PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);    
      if(partitionIDSource != partitionIDTarget){
	no_of_boundary_nodes++;
	block_bnd[partitionIDSource]++;   
	break; 
      }
    } endfor 
 } endfor
 return no_of_boundary_nodes;
}


double extended_quality_metrics::balance(graph_access& G) {
        std::vector<PartitionID> part_weights(G.get_partition_count(), 0);

        double overallWeight = 0;

        forall_nodes(G, n) {
                PartitionID curPartition = G.getPartitionIndex(n);
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

double extended_quality_metrics::balances(graph_access& G, std::vector<double> & block_balance) {
  std::vector<PartitionID> part_weights(G.get_partition_count(), 0);  
  double overallWeight = 0;
  forall_nodes(G, n){
    PartitionID curPartition = G.getPartitionIndex(n);
    (part_weights[curPartition])++;
    overallWeight++;
  } endfor
      
  double balance_part_weight = ceil(overallWeight / (double)G.get_partition_count());
  double cur_max             = -1;
  
  forall_blocks(G, p) {
    double cur = part_weights[p];
    if(cur > cur_max){
      cur_max = cur;
    }
    block_balance[p]=cur/balance_part_weight;
  } endfor

  double percentage = cur_max/balance_part_weight;
  return percentage;
} 

EdgeWeight extended_quality_metrics::objective(const PartitionConfig & config, graph_access & G, int* partition_map) {
        if(config.mh_optimize_communication_volume) {
                return max_communication_volume(G, partition_map);
        } else {
                return edge_cut(G, partition_map);
        }
}
