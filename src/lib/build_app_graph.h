#include <string>

#include "data_structure/graph_access.h"
#include "./karma_config.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "random_functions.h"
#include "timer.h"

int buildAppGraph(PartitionConfig & partition_config, std::string & graph_filename, graph_access & G);
