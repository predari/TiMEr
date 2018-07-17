#ifndef MAPPING_ALGORITHMS_H_
#define MAPPING_ALGORITHMS_H_

#include "configuration.h"
#include "data_structure/graph_access.h"
#include "partition/graph_partitioner.h"
#include "./../karma_config.h"
#include "./processor_graph_access.h"
#include "./alma_graph_access.h"
#include "./errors.h"

graph_access coarseComm (graph_access & graph, std::vector<int> & map, int parts);

graph_access removeLevel(graph_access & tree, unsigned int last);

void map(alma_graph_access & G, graph_access & C, processor_graph_access & P,
	 PartitionConfig & map_config, int64_t numberHierarchies, std::vector<int> & mapping);

void initialMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping);

void randomMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping);

void drbMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping);

void rcmMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping);

void multiSimpleMapping(graph_access & C, processor_graph_access & P, std::vector<int> & finalMapping);

void multiKahipMapping(graph_access & C, processor_graph_access & P, const std::string & type, std::vector<int> & finalMapping);

void brandfassMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping);

void brandfassLMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping);

void hoeflerMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping);

void hoeflerLMapping(graph_access & C, processor_graph_access & P, std::vector<int> & mapping);

void almaMapping(alma_graph_access &g, graph_access & C, processor_graph_access & P, PartitionConfig & map_config,
		 int64_t numberHierarchies, std::vector<int> & mapping);

#endif /* MAPPING_ALGORITHMS_H_ */
