#ifndef GRAPH_REPARTITIONER_H
#define GRAPH_REPARTITIONER_H

#include "graph_partitioner_m.h"

class graph_repartitioner : public graph_partitioner_m {
public:
private:
  virtual void single_run( PartitionConfig & config, graph_access & G);
};

#endif //include-guard
