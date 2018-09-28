/******************************************************************************
 * quality_metrics.h 
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

#include "data_structure/graph_access.h"
#include "./../karma_config.h"

class extended_quality_metrics {
public:
        extended_quality_metrics();
        virtual ~extended_quality_metrics ();

        EdgeWeight edge_cut(graph_access & G);
        EdgeWeight edge_cut(graph_access & G, int * partition_map); 
        EdgeWeight edge_cut(graph_access & G, PartitionID lhs, PartitionID rhs);
        EdgeWeight edgeCuts(graph_access & G, std::vector<EdgeWeight> & block_internal_edges, std::vector<EdgeWeight> & block_cut);
	void invade(graph_access & G, NodeID n, PartitionID pID, std::vector<bool> & visited);
	NodeID blockComps(graph_access & G);
        EdgeWeight max_communication_volume(graph_access & G, int * partition_map);
        EdgeWeight max_communication_volume(graph_access & G);
        EdgeWeight communicationVolumes(graph_access & G, std::vector<EdgeWeight> & block_volume);
        EdgeWeight objective(const PartitionConfig & config, graph_access & G, int * partition_map);
        int boundary_nodes(graph_access & G);
        int boundaryNodes(graph_access & G, std::vector<NodeID> & block_nodes, std::vector<NodeID> & block_bnd);
        double balance(graph_access & G);
	double balances(graph_access& G, std::vector<double> & block_balance);
};
