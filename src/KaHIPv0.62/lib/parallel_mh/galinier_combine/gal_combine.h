/******************************************************************************
 * gal_combine.h 
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

#ifndef GAL_COMBINE_XDMU5YB7
#define GAL_COMBINE_XDMU5YB7

#include "partition_config.h"
#include "data_structure/graph_access.h"

class gal_combine {
public:
        gal_combine();
        virtual ~gal_combine();

        void perform_gal_combine( PartitionConfig & config, graph_access & G);
};


#endif /* end of include guard: GAL_COMBINE_XDMU5YB7 */
