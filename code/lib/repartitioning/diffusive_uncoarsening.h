#ifndef DIFFUSIVE_UNCOARSENING_H
#define DIFFUSIVE_UNCOARSENING_H

#include "data_structure/graph_hierarchy.h"
#include "partition_config.h"

class diffusive_uncoarsening {
public:
  diffusive_uncoarsening( );
  virtual ~diffusive_uncoarsening();
  
  int perform_uncoarsening(const PartitionConfig & config, graph_hierarchy & hierarchy);
};


#endif //include-guard DIFFUSIVE_UNCOARSENING_H
