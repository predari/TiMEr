#include <cstdio>
#include <cmath>
#include <cstring> 
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <list>
#include <regex.h>
#include <sstream>
#include <string>

#include <argtable2.h>

#include "data_structure/graph_access.h"
#include "quality_metrics.h"
#include "graph_repartitioner.h"
#include "repartition_config.h"
#include "configuration.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "repartitioning/graph_partitioner.h"
#include "karma_config.h"
#include "random_functions.h"
#include "timer.h"
#include "repartitioning/cycle_improvements/cycle_refinement.h"

const int FAST           = 0;
const int ECO            = 1;
const int STRONG         = 2;
const int FASTSOCIAL     = 3;
const int ECOSOCIAL      = 4;
const int STRONGSOCIAL   = 5;

void prepPartitioning(PartitionConfig & partition_config, graph_access & G, bool repartition) {
  //code partially copy&pasted from kaffpa.cpp
       
  G.set_partition_count(partition_config.k); 
 
  NodeWeight largest_graph_weight = 0;
  forall_nodes(G, node) {
    largest_graph_weight += G.getNodeWeight(node);
  } 
  endfor
    
  double epsilon = (partition_config.imbalance)/100.0;
  if(  partition_config.imbalance == 0 ) {
    partition_config.upper_bound_partition      = (1+epsilon+0.01)*ceil(largest_graph_weight/(double)partition_config.k);
    partition_config.kaffpa_perfectly_balance   = true;
  } else {
    partition_config.upper_bound_partition      = (1+epsilon)*ceil(largest_graph_weight/(double)partition_config.k);
  }
  partition_config.largest_graph_weight       = largest_graph_weight;
  partition_config.graph_allready_partitioned = repartition;
  partition_config.kway_adaptive_limits_beta  = log(largest_graph_weight);

  std::cout <<  "block weight upper bound " <<  partition_config.upper_bound_partition  << std::endl;
        
  if(partition_config.input_partition != "") {
    std::cout <<  "reading input partition" << std::endl;
    graph_io::readPartition(G, partition_config.input_partition);
    partition_config.graph_allready_partitioned  = true;
    partition_config.no_new_initial_partitioning = true;
  }

  srand(partition_config.seed);
  random_functions::setSeed(partition_config.seed);

  std::cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << std::endl;
}

void partition(PartitionConfig& partition_config, graph_access& G, int k, int imbalance, int mode) {
  timer t;
  partition_config.k = k;
  configuration cfg;

  if (mode == FAST) cfg.fast(partition_config);
  else if (mode == ECO) cfg.eco(partition_config);
  else if (mode == STRONG) cfg.strong(partition_config);
  else if (mode == FASTSOCIAL) cfg.fastsocial(partition_config);
  else if (mode == ECOSOCIAL) cfg.ecosocial(partition_config);
  else if (mode == STRONGSOCIAL) cfg.strongsocial(partition_config);
  else {
    fprintf(stderr, "Invalid mode\n");
    exit(0);
  }

  partition_config.imbalance = imbalance;
  prepPartitioning(partition_config, G, false);
  graph_partitioner partitioner;
  quality_metrics qm;

  t.restart();
  partitioner.perform_partitioning(partition_config, G);
  std::cout <<  "time spent for partitioning " << t.elapsed()  << std::endl;
       
  // output some information about the partition that we have computed 
  std::cout << "cut \t\t"         << qm.edge_cut(G)                 << std::endl;
  std::cout << "finalobjective  " << qm.edge_cut(G)                 << std::endl;
  std::cout << "bnd \t\t"         << qm.boundary_nodes(G)           << std::endl;
  std::cout << "balance \t"       << qm.balance(G)                  << std::endl;
  std::cout << "finalbalance \t"  << qm.balance(G)                  << std::endl;
  std::cout << "max_comm_vol \t"  << qm.max_communication_volume(G) << std::endl;
}

void generateInput(PartitionConfig& partition_config, graph_access& G, int k, int imbalance, const char* filePath, int mode) {
  std::cout << "--- generate Input\n";
  timer t;
  graph_io::readGraphWeighted(G, filePath);
  std::cout << "io time: " << t.elapsed()  << std::endl;
  
  partition(partition_config, G, k, imbalance, mode);
}


void repartition(graph_access& G, int k, int imbalance, int mode) {
  std::cout << "--- repartition\n";

  RepartitionConfig repartition_config;

  repartition_config.k = k;
  configuration cfg;

  if (mode == FAST) cfg.fast(repartition_config);
  else if (mode == ECO) cfg.eco(repartition_config);
  else if (mode == STRONG) cfg.strong(repartition_config);
  else if (mode == FASTSOCIAL) cfg.fastsocial(repartition_config);
  else if (mode == ECOSOCIAL) cfg.ecosocial(repartition_config);
  else if (mode == STRONGSOCIAL) cfg.strongsocial(repartition_config);
  else {
    fprintf(stderr, "Invalid mode\n");
    exit(0);
  }

  repartition_config.use_fullmultigrid = false; //TODO: is this ok?
  repartition_config.use_wcycles = false;       //TODO: is this ok?

  repartition_config.label_propagation_refinement = false;
  repartition_config.quotient_graph_refinement_disabled = true;
  repartition_config.corner_refinement_enabled = false;
  repartition_config.softrebalance = false;
  repartition_config.refinement_type = REFINEMENT_TYPE_FM;//REFINEMENT_TYPE_FM_FLOW;

  repartition_config.imbalance = imbalance;
  
  repartition_config.diffusion_based_repartitioning = true;
  repartition_config.trunccons = true;
  repartition_config.trunccons_limited = true;
  repartition_config.trunccons_numDiffIters = 14;
  repartition_config.trunccons_numConsols = 9;
  repartition_config.trunccons_limited_useAllNeighbors = true;
  repartition_config.trunccons_limited_sizeOfNeighborhood = 8;

  repartition_config.scale_balancing = true;
  repartition_config.scale_balancing_iters = 16;
  repartition_config.flowBased_balancing = false;
  repartition_config.greedy_balancing = true;

  prepPartitioning(repartition_config, G, false);
  graph_repartitioner repartitioner;
  quality_metrics qm;

  timer t;
  t.restart();
  repartitioner.perform_partitioning(repartition_config, G);
  /*complete_boundary boundary(&G);
  boundary.build();
  cycle_refinement cr;
  cr.perform_refinement(partition_config, G, boundary);*/
  std::cout <<  "time spent for repartitioning " << t.elapsed()  << std::endl;
       
  // output some information about the partition that we have computed 
  std::cout << "cut \t\t"         << qm.edge_cut(G)                 << std::endl;
  std::cout << "finalobjective  " << qm.edge_cut(G)                 << std::endl;
  std::cout << "bnd \t\t"         << qm.boundary_nodes(G)           << std::endl;
  std::cout << "balance \t"       << qm.balance(G)                  << std::endl;
  std::cout << "finalbalance \t"  << qm.balance(G)                  << std::endl;
  std::cout << "max_comm_vol \t"  << qm.max_communication_volume(G) << std::endl;
}

void fromUnbalanced(const char* file, int k, int genImbalance, int reImbalance, int mode, const char*  file_out) {
  std::cout << "File: \t\t\t\t" << file << std::endl;
  std::cout << "#Partitions: \t\t\t" << k << std::endl;
  std::cout << "generator Imbalance: \t\t" << genImbalance << std::endl;
  std::cout << "repartitioner Imbalance: \t" << reImbalance << std::endl;
  std::cout << "mode: \t\t\t\t" << mode << std::endl;

  PartitionConfig in_config;
  graph_access G;
  generateInput(in_config, G, k, genImbalance, file, mode);

  repartition(G, k, reImbalance, mode);

  if (file_out != NULL) graph_io::writePartition(G, file_out);
}

void fromSeqGenerator(const char* prefix, int N, int k, int imbalance, int mode, const char* file_out) {
  std::cout << "--- Import data from sequence generator\n";
  
  //read and import first graph
  FILE* f;
  char* fileG1 = new char[strlen(prefix)+12];
  strcpy(fileG1, prefix);
  strncat(fileG1, "-00000.graph", 12);
  graph_access* G = new graph_access();
  
  graph_io::readGraphWeighted(*G, fileG1);

  //partition first graph

  std::cout << "--- partition initial graph\n";
  PartitionConfig G_config;
  partition(G_config, *G, k, imbalance, mode);

  //read and import next frames

  for (int frame = 1; frame <= N; ++frame) {
    graph_access* G_tmp = new graph_access();
    std::stringstream sst;
    sst << std::setfill('0') << std::setw(5) << frame;
    std::string gtfile = std::string(prefix) + "-" + sst.str() + ".graph";
    std::string atfile = std::string(prefix) + "-" + sst.str() + ".at";

    graph_io::readGraphWeighted(*G_tmp, gtfile.c_str());

    //read ancestors and transform partitioning
    int* anc;
    f = fopen(atfile.c_str(), "r");
    anc = new int[G_tmp->number_of_nodes()];
    for (int i = 0; i < G_tmp->number_of_nodes(); ++i) {
      fscanf(f, "%d", &anc[i]);
    }
    fclose(f);

    forall_nodes((*G_tmp), n) {
      G_tmp->setPartitionIndex(n, G->getPartitionIndex(anc[n]));
    } endfor
        
    delete G;
    G = G_tmp;
  }

  quality_metrics qm;
  std::cout << "--- New Graph with old partition:\n";
  std::cout << "cut \t\t"         << qm.edge_cut(*G)                 << std::endl;
  std::cout << "bnd \t\t"         << qm.boundary_nodes(*G)           << std::endl;
  std::cout << "balance \t"       << qm.balance(*G)                  << std::endl;
  std::cout << "max_comm_vol \t"  << qm.max_communication_volume(*G) << std::endl;
      

  //repartition
  repartition(*G, k, imbalance, mode);
      
  if (file_out != NULL) graph_io::writePartition(*G, file_out);
}

int main(int argn, char** argv) {
  const char *progname = argv[0];
  
  // Setup argtable parameters.
  
  struct arg_rex *g_cmd = arg_rex1(NULL, NULL, "gen", NULL, REG_ICASE, NULL);
  struct arg_str *g_filename = arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file to partition.");
  struct arg_dbl *g_genImbalance = arg_dbl0(NULL, "genImbalance", NULL, "Desired balance. Default: 20 (%).");

  struct arg_rex *s_cmd = arg_rex1(NULL, NULL, "seq", NULL, REG_ICASE, NULL);
  struct arg_str *s_prefix = arg_strn(NULL, NULL, "PREFIX", 1, 1, "Prefix to sequence files.");
  struct arg_int *s_N = arg_int0(NULL, "N", NULL, "Number of frames.");

  struct arg_str *sg_filename_output = arg_str0(NULL, "output_filename", NULL, "Specify the name of the output file (that contains the partition).");
  struct arg_dbl *sg_imbalance = arg_dbl0(NULL, "imbalance", NULL, "Desired balance. Default: 3 (%).");
  struct arg_rex *sg_preconfiguration = arg_rex0(NULL, "preconfiguration", "^(strong|eco|fast|fastsocial|ecosocial|strongsocial)$", "VARIANT", REG_EXTENDED, "Use a preconfiguration. (Default: strong) [strong|eco|fast|fastsocial|ecosocial|strongsocial]." );
  struct arg_lit *sg_help = arg_lit0(NULL, "help","Print help.");
  struct arg_int *sg_k = arg_int1(NULL, "k", NULL, "Number of blocks to partition the graph.");

  struct arg_end *end = arg_end(100);

  // Define argtable.
  void* g_argtable[] = {g_cmd, sg_help, g_filename, g_genImbalance, sg_filename_output, sg_imbalance, sg_preconfiguration, sg_k, end};
  void* s_argtable[] = {s_cmd, sg_help, s_prefix, s_N, sg_filename_output, sg_imbalance, sg_preconfiguration, sg_k, end};

  // Parse arguments.
  int g_nerrors = arg_parse(argn, argv, g_argtable);
  int s_nerrors = arg_parse(argn, argv, s_argtable);

  bool seq = true;
  bool gen = true;

  if (g_nerrors > 0) {
    gen = false; 
  }

  if (s_nerrors > 0) {
    seq = false;
  }

  if (!gen and !seq) {
    arg_print_errors(stderr, end, progname);
    printf("Try '%s --help' for more information.\n",progname);
    arg_freetable(g_argtable, sizeof(g_argtable) / sizeof(g_argtable[0]));
    arg_freetable(s_argtable, sizeof(s_argtable) / sizeof(s_argtable[0]));
    return 1; 
  }

  // Catch case that help was requested.
  if (sg_help->count > 0) {
    printf("Usage: %s", progname);
    //TODO: test this
    if (seq) {
      arg_print_syntax(stdout, s_argtable, "\n");
      arg_print_glossary(stdout, s_argtable,"  %-40s %s\n");
    }
    if (gen) {
      arg_print_syntax(stdout, g_argtable, "\n");
      arg_print_glossary(stdout, g_argtable,"  %-40s %s\n");
    }
    arg_freetable(g_argtable, sizeof(g_argtable) / sizeof(g_argtable[0]));
    arg_freetable(s_argtable, sizeof(s_argtable) / sizeof(s_argtable[0]));
    return 1;
  }
  
  int k = 16;
  int imbalance = 3;
  int mode = STRONG;
  const char* file_out = "tmppartition";

  if (sg_k->count > 0) {
    k = sg_k->ival[0];
  }

  if(sg_preconfiguration->count > 0) {
    if(strcmp("strong", sg_preconfiguration->sval[0]) == 0) {
      mode = STRONG;
    } else if (strcmp("eco", sg_preconfiguration->sval[0]) == 0) {
      mode = ECO;
    } else if (strcmp("fast", sg_preconfiguration->sval[0]) == 0) {
      mode = FAST;
    } else if (strcmp("fastsocial", sg_preconfiguration->sval[0]) == 0) {
      mode = FASTSOCIAL;
    } else if (strcmp("ecosocial", sg_preconfiguration->sval[0]) == 0) {
      mode = ECOSOCIAL;
    } else if (strcmp("strongsocial", sg_preconfiguration->sval[0]) == 0) {
      mode = STRONGSOCIAL;
    } else {
      fprintf(stderr, "Invalid preconfiguration variant: \"%s\"\n", sg_preconfiguration->sval[0]);
      exit(0);
    }
  }

  if(sg_filename_output->count > 0) {
    file_out = sg_filename_output->sval[0];
  }

  if(sg_imbalance->count > 0) {
    imbalance = sg_imbalance->dval[0];
  }

  if (seq) {
    const char* prefix = "../mesh-sequences/refine-";
    int N = 100;
    if(s_prefix->count > 0) {
      prefix = s_prefix->sval[0];
    }
    if(s_N->count > 0) {
      N = s_N->ival[0];
    }
    fromSeqGenerator(prefix, N, k, imbalance, mode, file_out);
  }

  if (gen) {
    int genImbalance = 20;
    const char* file = "../graph_part_viz/samplegraphs/Grid100x100";

    if(g_genImbalance->count > 0) {
      genImbalance = g_genImbalance->dval[0];
    }

    if(g_filename->count > 0) {
      file = g_filename->sval[0];
    }
    fromUnbalanced(file, k, genImbalance, imbalance, mode, file_out);
  }

  return 0;
}
