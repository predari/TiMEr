/******************************************************************************
 * graph_partitioner_m.h 
 *
 * copied and slightely modified from graph_partitioner.h:
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

#ifndef GRAPH_PARTITIONER_M_H
#define GRAPH_PARTITIONER_M_H

#include "coarsening/coarsening.h"
#include "coarsening/stop_rules/stop_rules.h"
#include "data_structure/graph_access.h"
#include "../karma_config.h"
#include "uncoarsening/refinement/refinement.h"

class graph_partitioner_m {
public:
        graph_partitioner_m();
        virtual ~graph_partitioner_m();

        virtual void perform_partitioning(PartitionConfig & graph_partitioner_m_config, graph_access & G);
        virtual void perform_recursive_partitioning(PartitionConfig & graph_partitioner_m_config, graph_access & G);

private:
        virtual void perform_recursive_partitioning_internal(PartitionConfig & graph_partitioner_m_config, 
                                                     graph_access & G, 
                                                     PartitionID lb, PartitionID ub);
        virtual void single_run( PartitionConfig & config, graph_access & G);

        unsigned m_global_k;
	int m_global_upper_bound;
        int m_rnd_bal;
};

#endif /* end of include guard: GRAPH_PARTITIONER_M_H */
