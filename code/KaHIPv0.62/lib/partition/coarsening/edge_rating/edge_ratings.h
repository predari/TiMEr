/******************************************************************************
 * edge_ratings.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef EDGE_RATING_FUNCTIONS_FUCW7H6Y
#define EDGE_RATING_FUNCTIONS_FUCW7H6Y

#define SQR_NR_BINS 10
#define MAX_NUMBER_TREES 20
#define ONE_OVER_UPPER_QUANTILE 0//20 means upper 5 percent 

#include "data_structure/graph_access.h"
#include "partition_config.h"       

typedef struct DirEdge{
  NodeID source;
  NodeID target;
  EdgeID edID;
  struct DirEdge* next;
} DirEdge;

class edge_ratings {
public:
        edge_ratings(const PartitionConfig & partition_config);
        virtual ~edge_ratings();

        void rate(graph_access & G, unsigned level);
        void rate_weight(graph_access & G);
        void rate_expansion_star_2(graph_access & G);
        void rate_expansion_star(graph_access & G);
        void rate_expansion_star_2_algdist(graph_access & G);
        void rate_inner_outer(graph_access & G);
        void rate_pseudogeom(graph_access & G); 
        void compute_algdist(graph_access & G, std::vector<float> & dist); 
        void rate_zero(graph_access & G);

	/******************************* djoko ******************************/
	DirEdge* buildSortedEdgeList(graph_access & G, std::vector<unsigned int> & bottleneckIndex1, unsigned int maxBI1);
	void dfs(NodeID n, graph_access & G, std::vector<bool> & marked);
	NodeID getNumberConnectedComponents(graph_access & G);
	double getLocalClusteringCoefficient(graph_access & G);
	void characterizeG(graph_access & G, bool cumul);
	double rand01(void);
	EdgeID getEdge(NodeID s, NodeID t, graph_access & G);
	void getReverseEdges(graph_access & G, std::vector<EdgeID> & revEdge);
	void globalBFS(NodeID s, graph_access & G, std::vector<bool> & inBFSTree);
	void parentMaker(NodeID w, graph_access & G, std::vector<bool> & inMST);
	void globalMST(NodeID MSTroot, graph_access & G, std::vector<EdgeID> & revEdge, std::vector<unsigned int> & bottleneckIndex1, unsigned int maxBI1, 
		       std::vector<bool> & inMST);
	NodeID commonAncestor(NodeID v1, NodeID v2,
			      std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants);
	void resetInTree(graph_access & G, std::vector<bool> & inTree);
	void getLabels(NodeID w, graph_access & G,
		       std::vector<bool> & inMST, NodeID* currLabel_ptr,
		       std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants);
	void leaf(NodeID v, graph_access & G, std::vector<bool> & inMST,
		  std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants,
		  std::vector<double> & intraWeight, std::vector<double> & interWeight,
		  std::vector<double> & subtreeWeight);
	void update(NodeID u, graph_access & G, std::vector<bool> & inMST, 
		    std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants,
		    std::vector<double> & intraWeight, std::vector<double> & interWeight,
		    std::vector<double> & subtreeWeight);
	void mainDFT(NodeID v, graph_access & G, std::vector<bool> & inMST, 
		     std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants,
		     std::vector<double> & intraWeight, std::vector<double> & interWeight,
		     std::vector<double> & subtreeWeight);
	unsigned int computeBottleneckIndices_one(graph_access & G,
					    std::vector<EdgeID> & revEdge,
						  std::vector<unsigned int> & bottleneckIndex1);
	EdgeID findProxyEdge1(graph_access & G,
			      std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants,
			      std::vector<bool> & inMST, std::vector<unsigned int> & bottleneckIndex1,
			      NodeID source, NodeID target);
	EdgeID findProxyEdgeR(graph_access & G,
			      std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants,
			      std::vector<bool> & inMST, std::vector<EdgeRatingType> & bottleneckIndexR,
			      NodeID source, NodeID target);
  
	void computeBottleneckIndices_resistance(graph_access & G,
						 std::vector<EdgeID> & revEdge,
						 std::vector<unsigned int> & bottleneckIndex1,
						 unsigned int maxBI1,
						 std::vector<EdgeRatingType> & bottleneckIndexR);
	void writeBottleneckIndexHistoR(graph_access & G, std::vector<EdgeRatingType> & bottleneckIndexR, bool cumul);
	void rate_resistance(graph_access & G, unsigned level); 
	/********************************************************************/
	
 private:
        const PartitionConfig & partition_config;
};


#endif /* end of include guard: EDGE_RATING_FUNCTIONS_FUCW7H6Y */
