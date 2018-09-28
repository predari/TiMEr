/******************************************************************************
 * parse_parameters.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013 Christian Schulz <christian.schulz@kit.edu>
 * Original version from Christian Schulz <christian.schulz@kit.edu>
 *
 * Copied and adapted by Roland Glantz <roland.glantz@kit.edu>
 * All changes are additions (see "parse_repart_parameters()" and the blocks
 * that apply to the "MODE_PART", "MODE_REPART" and "MODE_MAP" directives
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


#ifndef PARSE_PARAMETERS_GPJMGSM8
#define PARSE_PARAMETERS_GPJMGSM8

#include <omp.h>
#include "./../lib/karma_config.h"
#include "configuration.h"

int parse_parameters(int argn, char **argv,
                     PartitionConfig & partition_config,
                     std::string & graph_filename, 
                     bool & is_graph_weighted, 
                     bool & suppress_program_output, 
                     bool & recursive) {

  const char *progname = argv[0];

  // Setup argtable parameters.
  struct arg_lit *help                                 = arg_lit0(NULL, "help","Print help.");
  struct arg_lit *edge_rating_tiebreaking              = arg_lit0(NULL, "edge_rating_tiebreaking","Enable random edgerating tiebreaking.");
  struct arg_lit *match_islands                        = arg_lit0(NULL, "match_islands","Enable matching of islands during gpa algorithm.");
  struct arg_lit *only_first_level                     = arg_lit0(NULL, "only_first_level","Disable Multilevel Approach. Only perform on the first level. (Currently only initial partitioning).");
  struct arg_lit *graph_weighted                       = arg_lit0(NULL, "weighted","Read the graph as weighted graph.");
  struct arg_lit *enable_corner_refinement             = arg_lit0(NULL, "enable_corner_refinement","Enables corner refinement.");
  struct arg_lit *disable_qgraph_refinement            = arg_lit0(NULL, "disable_qgraph_refinement","Disables qgraph refinement.");
  struct arg_lit *use_fullmultigrid                    = arg_lit0(NULL, "use_fullmultigrid","Enable full multigrid (wcycles have to be enabled).");
  struct arg_lit *use_vcycle                           = arg_lit0(NULL, "use_vcycle","Enable vcycle .");
  struct arg_lit *compute_vertex_separator             = arg_lit0(NULL, "compute_vertex_separator","Compute vertex separator.");
  struct arg_lit *first_level_random_matching          = arg_lit0(NULL, "first_level_random_matching", "The first level will be matched randomly.");
  struct arg_lit *rate_first_level_inner_outer         = arg_lit0(NULL, "rate_first_level_inner_outer", "The edge rating for the first level is inner outer.");
  struct arg_lit *use_bucket_queues                    = arg_lit0(NULL, "use_bucket_queues", "Use bucket priority queues during refinement.");
  struct arg_lit *use_wcycles                          = arg_lit0(NULL, "use_wcycle", "Enables wcycles.");
  struct arg_lit *disable_refined_bubbling             = arg_lit0(NULL, "disable_refined_bubbling", "Disables refinement during initial partitioning using bubbling (Default: enabled).");
  struct arg_lit *enable_convergence                   = arg_lit0(NULL, "enable_convergence", "Enables convergence mode, i.e. every step is running until no change.(Default: disabled).");
  struct arg_lit *enable_omp                           = arg_lit0(NULL, "enable_omp", "Enable parallel omp.");
  struct arg_lit *wcycle_no_new_initial_partitioning   = arg_lit0(NULL, "wcycle_no_new_initial_partitioning", "Using this option, the graph is initially partitioned only the first time we are at the deepest level.");
  struct arg_str *filename                             = arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file to partition.");
  struct arg_str *filename_output                      = arg_str0(NULL, "filename_output", NULL, "Specify the name of the output file (that contains a partition or a mapping).");
  struct arg_int *user_seed                            = arg_int0(NULL, "seed", NULL, "Seed to use for the PRNG.");
  struct arg_int *k                                    = arg_int0(NULL, "k", NULL, "Number of blocks to partition the graph.");
  struct arg_rex *edge_rating                          = arg_rex0(NULL, "edge_rating", "^(weight|expansionstar|expansionstar2|expansionstar2algdist|geom|resistance|zero)$", "RATING", REG_EXTENDED, "Edge rating to use. One of {weight|expansionstar|expansionstar2|expansionstar2algdist|geom|resistance|zero}. Default: weight");
  struct arg_rex *refinement_type                      = arg_rex0(NULL, "refinement_type", "^(fm|fm_flow|flow)$", "TYPE", REG_EXTENDED, "Refinementvariant to use. One of {fm, fm_flow, flow}. Default: fm"  );
  struct arg_rex *matching_type                        = arg_rex0(NULL, "matching", "^(random|hem|shem|regions|gpa|randomgpa|localmax)$", "TYPE", REG_EXTENDED, "Type of matchings to use during coarsening. One of {random, hem," " shem, regions, gpa, randomgpa, localmax}."  );
  struct arg_int *mh_pool_size                         = arg_int0(NULL, "mh_pool_size", NULL, "MetaHeuristic Pool Size.");
  struct arg_lit *mh_plain_repetitions                 = arg_lit0(NULL, "mh_plain_repetitions", "");
  struct arg_lit *mh_penalty_for_unconnected           = arg_lit0(NULL, "mh_penalty_for_unconnected", "Add a penalty on the objective function if the computed partition contains blocks that are not connected.");
  struct arg_lit *mh_disable_nc_combine                = arg_lit0(NULL, "mh_disable_nc_combine", "");
  struct arg_lit *mh_disable_cross_combine             = arg_lit0(NULL, "mh_disable_cross_combine", "");
  struct arg_lit *mh_disable_combine                   = arg_lit0(NULL, "mh_disable_combine", "");
  struct arg_lit *mh_enable_quickstart                 = arg_lit0(NULL, "mh_enable_quickstart", "Enables the quickstart option.");
  struct arg_lit *mh_disable_diversify_islands         = arg_lit0(NULL, "mh_disable_diversify_islands", "");
  struct arg_lit *mh_disable_diversify                 = arg_lit0(NULL, "mh_disable_diversify", "");
  struct arg_lit *mh_diversify_best                    = arg_lit0(NULL, "mh_diversify_best", "Uses best individuum instead of random during diversification.");
  struct arg_lit *mh_enable_tournament_selection       = arg_lit0(NULL, "mh_enable_tournament_selection", "Enables the tournament selection roule instead of choosing two random inidiviuums.");
  struct arg_lit *mh_cross_combine_original_k          = arg_lit0(NULL, "mh_cross_combine_original_k", "");
  struct arg_lit *mh_optimize_communication_volume     = arg_lit0(NULL, "mh_optimize_communication_volume", "Fitness function is modified to optimize communication volume instead of the number of cut edges.");
  struct arg_lit *disable_balance_singletons           = arg_lit0(NULL, "disable_balance_singletons", "");
  struct arg_lit *gpa_grow_internal                    = arg_lit0(NULL, "gpa_grow_internal", "If the graph is allready partitions the paths are grown only block internally.");
  struct arg_int *initial_partitioning_repetitions     = arg_int0(NULL, "initial_partitioning_repetitions", NULL, "Number of initial partitioning repetitons. Default: 5.");
  struct arg_int *minipreps                            = arg_int0(NULL, "minipreps", NULL, "Default: 10.");
  struct arg_int *aggressive_random_levels             = arg_int0(NULL, "aggressive_random_levels", NULL, "In case matching is randomgpa, this is the number of levels that should be matched using random matching. Default: 3.");
  struct arg_dbl *imbalance                            = arg_dbl0(NULL, "imbalance", NULL, "Desired balance. Default: 3 (%).");
  struct arg_rex *initial_partition                    = arg_rex0(NULL, "initial_partitioner", "^(metis|scotch|hybrid|bubbling|squeez|metaheuristic|recursive)$", "PARTITIONER", REG_EXTENDED, "Type of matchings to use during coarsening. One of {metis, scotch, bubbling, hybrid, recursive)." );
  struct arg_lit *initial_partition_optimize           = arg_lit0(NULL, "initial_partition_optimize", "Enables postoptimization of initial partition.");
  struct arg_rex *bipartition_algorithm                = arg_rex0(NULL, "bipartition_algorithm", "^(bfs|fm|squeezing)$", "TYPE", REG_EXTENDED, "Type of bipartition algorithm to use in case of recursive partitioning. One of " " {bfs, fm, squeezing}."  );
  struct arg_rex *permutation_quality                  = arg_rex0(NULL, "permutation_quality", "^(none|fast|good|cacheefficient)$", "QUALITY", REG_EXTENDED, "The quality of permutations to use. One of {none, fast," " good, cacheefficient}."  );
  struct arg_rex *permutation_during_refinement        = arg_rex0(NULL, "permutation_during_refinement", "^(none|fast|good)$", "QUALITY", REG_EXTENDED, "The quality of permutations to use during 2way fm refinement. One of {none, fast," " good}."  );
  struct arg_int *fm_search_limit                      = arg_int0(NULL, "fm_search_limit", NULL, "Search limit for 2way fm local search: Default 1 (%).");
  struct arg_int *bipartition_post_fm_limit            = arg_int0(NULL, "bipartition_post_fm_limit", NULL, "Search limit for the fm search after a bipartition has been created. :");
  struct arg_int *bipartition_post_ml_limit            = arg_int0(NULL, "bipartition_post_ml_limit", NULL, "Search limit for the multilevel fm search after a bipartition has been created. :");
  struct arg_int *bipartition_tries                    = arg_int0(NULL, "bipartition_tries", NULL, "Number of tries to find a bipartition (during recursive intial partitioning).");
  struct arg_rex *refinement_scheduling_algorithm      = arg_rex0(NULL, "refinement_scheduling_algorithm", "^(fast|active_blocks|active_blocks_kway)$", "QUALITY", REG_EXTENDED, " One of {fast, active_blocks, active_blocks_kway}.");
  struct arg_dbl *bank_account_factor                  = arg_dbl0(NULL, "bank_account_factor", NULL, "The bank account factor for the scheduler. Default 1.5 (%).");
  struct arg_dbl *flow_region_factor                   = arg_dbl0(NULL, "flow_region_factor", NULL, "If using flow, then the regions found are sized flow_region_factor * imbalance. Default: 4 (%).");
  struct arg_dbl *kway_adaptive_limits_alpha           = arg_dbl0(NULL, "kway_adaptive_limits_alpha", NULL, "This is the factor alpha used for the adaptive stopping criteria. Default: 1.0");
  struct arg_rex *stop_rule                            = arg_rex0(NULL, "stop_rule", "^(simple|multiplek|strong)$", "VARIANT", REG_EXTENDED, "Stop rule to use. One of {simple, multiplek, strong}. Default: simple" );
  struct arg_int *num_vert_stop_factor                 = arg_int0(NULL, "num_vert_stop_factor", NULL, "x*k (for multiple_k stop rule). Default 20.");
  struct arg_rex *kway_search_stop_rule                = arg_rex0(NULL, "kway_stop_rule", "^(simple|adaptive)$", "VARIANT", REG_EXTENDED, "Stop rule to use during kway_refinement. One of {simple, adaptive}. Default: simple" );
  struct arg_int *bubbling_iterations                  = arg_int0(NULL, "bubbling_iterations", NULL, "Number of bubbling iterations to perform: Default 1 .");
  struct arg_int *kway_rounds                          = arg_int0(NULL, "kway_rounds", NULL, "Number of kway refinement rounds to perform: Default 1 .");
  struct arg_int *kway_fm_limits                       = arg_int0(NULL, "kway_fm_search_limit", NULL, "Search limit for kway fm local search: Default 1 .");
  struct arg_int *global_cycle_iterations              = arg_int0(NULL, "global_cycle_iterations", NULL, "Number of V-cycle iterations: Default 2.");
  struct arg_int *level_split                          = arg_int0(NULL, "level_split", NULL, "Number of trial tree levels (1 means on each level a two trials are performed). Default: 2.");
  struct arg_int *toposort_iterations                  = arg_int0(NULL, "toposort_iterations", NULL, "Number of topo sort iterations). Default: 4.");
  struct arg_lit *most_balanced_flows                  = arg_lit0(NULL, "most_balanced_flows", "(Default: disabled)");
  struct arg_str *input_partition                      = arg_str0(NULL, "input_partition", NULL, "Input partition to use.");
  struct arg_lit *recursive_bipartitioning             = arg_lit0(NULL, "recursive_bipartitioning", "Use recursive bipartitioning instead of kway methods.");
  struct arg_lit *suppress_output                      = arg_lit0(NULL, "suppress_output", "(Default: output enabled)");
  struct arg_lit *disable_max_vertex_weight_constraint = arg_lit0(NULL, "disable_max_vertex_weight_constraint", "Disables the max vertex weight constraint during the contraction.");
  struct arg_int *local_multitry_fm_alpha              = arg_int0(NULL, "local_multitry_fm_alpha", NULL, "Search limit factor alpha for multitry fm.");
  struct arg_int *local_multitry_rounds                = arg_int0(NULL, "local_multitry_rounds", NULL, "Number of rounds for local multitry fm.");
  struct arg_int *initial_partition_optimize_fm_limits = arg_int0(NULL, "initial_partition_optimize_fm_limits", NULL, "Initial Partition Optimize FM limits. (Default: 20)");
  struct arg_int *initial_partition_optimize_multitry_fm_alpha = arg_int0(NULL, "initial_partition_optimize_multitry_fm_limits", NULL, "Initial Partition Optimize Multitry FM limits. (Default: 20)");
  struct arg_int *initial_partition_optimize_multitry_rounds   = arg_int0(NULL, "initial_partition_optimize_multitry_rounds", NULL, "(Default: 100)");

#ifdef MODE_KAFFPA
  struct arg_rex *preconfiguration                     = arg_rex0(NULL, "preconfiguration", "^(strong|eco|fast|fastsocial|ecosocial|strongsocial)$", "VARIANT", REG_EXTENDED, "Use a preconfiguration. (Default: eco) [strong|eco|fast|fastsocial|ecosocial|strongsocial]." );
#endif

#ifdef MODE_GRAPH_INFO
  struct arg_rex *preconfiguration                     = arg_rex0(NULL, "preconfiguration", "^(strong|eco|fast|fastsocial|ecosocial|strongsocial)$", "VARIANT", REG_EXTENDED, "Use a preconfiguration. (Default: eco) [strong|eco|fast|fastsocial|ecosocial|strongsocial]." );
#endif

#ifdef MODE_PARTITION_INFO
  struct arg_rex *preconfiguration                     = arg_rex0(NULL, "preconfiguration", "^(strong|eco|fast|fastsocial|ecosocial|strongsocial)$", "VARIANT", REG_EXTENDED, "Use a preconfiguration. (Default: eco) [strong|eco|fast|fastsocial|ecosocial|strongsocial]." );
#endif

#ifdef MODE_PART
  struct arg_rex *preconfiguration                     = arg_rex0(NULL, "preconfiguration", "^(strong|eco|fast|fastsocial|ecosocial|strongsocial)$", "VARIANT", REG_EXTENDED, "Use a preconfiguration. (Default: eco) [strong|eco|fast|fastsocial|ecosocial|strongsocial]." );
#endif

#ifdef MODE_REPART
  struct arg_rex *preconfiguration                     = arg_rex0(NULL, "preconfiguration", "^(strong|eco|fast|fastsocial|ecosocial|strongsocial)$", "VARIANT", REG_EXTENDED, "Use a preconfiguration. (Default: eco) [strong|eco|fast|fastsocial|ecosocial|strongsocial]." );
#endif

#ifdef MODE_MAP
  struct arg_rex *preconfiguration                     = arg_rex0(NULL, "preconfiguration", "^(strong|eco|fast|fastsocial|ecosocial|strongsocial)$", "VARIANT", REG_EXTENDED, "Use a preconfiguration. (Default: eco) [strong|eco|fast|fastsocial|ecosocial|strongsocial]." );

  struct arg_rex *multikahip_preconfiguration          = arg_rex0(NULL, "multikahip_preconfiguration", "^(strong|eco|fast|fastsocial|ecosocial|strongsocial)$", "VARIANT", REG_EXTENDED, "[strong|eco|fast|fastsocial|ecosocial|strongsocial] (default eco)");
#endif

#ifdef MODE_KAFFPAE
  struct arg_rex *preconfiguration                     = arg_rex0(NULL, "preconfiguration", "^(strong|eco|fast|fastsocial|ecosocial|strongsocial)$", "VARIANT", REG_EXTENDED, "Use a preconfiguration. (Default: strong) [strong|eco|fast|fastsocial|ecosocial|strongsocial]." );
#endif

#ifdef MODE_PARTITIONTOVERTEXSEPARATOR
  struct arg_rex *preconfiguration                     = arg_rex0(NULL, "preconfiguration", "^(strong|eco|fast|fastsocial|ecosocial|strongsocial)$", "VARIANT", REG_EXTENDED, "Use a preconfiguration. (Default: strong) [strong|eco|fast|fastsocial|ecosocial|strongsocial]." );
#endif

#ifdef MODE_LABELPROPAGATION
  struct arg_rex *preconfiguration                     = arg_rex0(NULL, "preconfiguration", "^(strong|eco|fast|fastsocial|ecosocial|strongsocial)$", "VARIANT", REG_EXTENDED, "Use a preconfiguration. (Default: strong) [strong|eco|fast|fastsocial|ecosocial|strongsocial]." );
#endif

  struct arg_dbl *time_limit                           = arg_dbl0(NULL, "time_limit", NULL, "Time limit in s. Default 0s .");
  struct arg_int *unsuccessful_reps                    = arg_int0(NULL, "unsuccessful_reps", NULL, "Unsuccessful reps to fresh start.");
  struct arg_int *local_partitioning_repetitions       = arg_int0(NULL, "local_partitioning_repetitions", NULL, "Number of local repetitions.");
  struct arg_int *amg_iterations                       = arg_int0(NULL, "amg_iterations", NULL, "Number of amg iterations.");
  struct arg_int *mh_flip_coin                         = arg_int0(NULL, "mh_flip_coin", NULL, "Control the ratio of mutation and crossovers. c/10 Mutation and (10-c)/10 crossovers.");
  struct arg_int *mh_initial_population_fraction       = arg_int0(NULL, "mh_initial_population_fraction", NULL, "Control the initial population fraction parameter (Default: 1000).");
  struct arg_lit *mh_print_log                         = arg_lit0(NULL, "mh_print_log", "Each PE prints a logfile (timestamp, edgecut).");
  struct arg_lit *mh_sequential_mode                   = arg_lit0(NULL, "mh_sequential_mode", "Disables all MH algorithms. Use KaFFPa in a parallel setting.");
  struct arg_rex *kaba_neg_cycle_algorithm             = arg_rex0(NULL, "kaba_neg_cycle_algorithm", "^(ultramodel|randomcycle|playfield|ultramodelplus)$", "VARIANT", REG_EXTENDED, "Balanced refinement operator to use. On of randomcycle, ultramodel, playfield, ultramodelplus" );
  struct arg_dbl *kabaE_internal_bal                   = arg_dbl0(NULL, "kabaE_internal_bal", NULL, "Control the internal balance paramter for kaffpaE (Default: 0.01) (1 percent)");
  struct arg_int *kaba_internal_no_aug_steps_aug       = arg_int0(NULL, "kaba_internal_no_aug_steps_aug", NULL, "Internal number of steps in the augmented models of negative cycle detection.");
  struct arg_int *kaba_packing_iterations              = arg_int0(NULL, "kaba_packing_iterations", NULL, "Number of packing iterations.");
  struct arg_int *kaba_unsucc_iterations               = arg_int0(NULL, "kaba_unsucc_iterations", NULL, "Number of unsucc iterations until a rebalancing step is performed.");
  struct arg_lit *kaba_flip_packings                   = arg_lit0(NULL, "kaba_flip_packings", "Enable flip packing mode (if ultramodelplus is used).");
  struct arg_rex *kaba_lsearch_p                       = arg_rex0(NULL, "kaba_lsearch_p", "^(coindiff|coinrnd|nocoindiff|nocoinrnd)$", "VARIANT", REG_EXTENDED, "Make more localized search in ultraplus model.");
  struct arg_lit *kaffpa_perfectly_balanced_refinement = arg_lit0(NULL, "kaffpa_perfectly_balanced_refinement", "Enable perfectly balanced refinement during ML KaFFPa.");
  struct arg_lit *kaba_disable_zero_weight_cycles      = arg_lit0(NULL, "kaba_disable_zero_weight_cycles", "Disable zero weight cycle diversification.");
  struct arg_lit *enforce_balance                      = arg_lit0(NULL, "enforce_balance", "Uses eps+1 to create a partition, and then runs rebalancing and negative cycle detection to output a partition that fulfills the eps-balance constraint.");
  struct arg_lit *mh_enable_tabu_search                = arg_lit0(NULL, "mh_enable_tabu_search", "Enables our version of combine operation by block matching; +tabusearch and all our refinement algorithms.");
  struct arg_lit *mh_enable_kabapE                     = arg_lit0(NULL, "mh_enable_kabapE", "Enable combine operator KaBaPE");
  struct arg_int *maxT                                 = arg_int0(NULL, "maxT", NULL, "maxT parameter for Tabu Search");
  struct arg_int *maxIter                              = arg_int0(NULL, "maxIter", NULL, "maxIter parameter for Tabu Search");

  struct arg_int *cluster_upperbound                   = arg_int0(NULL, "cluster_upperbound", NULL, "Set a size-constraint on the size of a cluster. Default: none");
  struct arg_int *label_propagation_iterations         = arg_int0(NULL, "label_propagation_iterations", NULL, "Set the number of label propgation iterations. Default: 10.");

#ifdef MODE_REPART
  struct arg_str *filename_old_graph                   = arg_str0(NULL, "filename_old_graph", NULL, "Path to file of old graph");
  struct arg_str *filename_ancestors                   = arg_str0(NULL, "filename_ancestors", NULL, "Ancestors of nodes in new graph");
  struct arg_dbl *new_imbalance                        = arg_dbl1(NULL,"new_imbalance", NULL,"upper bound for imbalance of repartition");
  struct arg_int *new_k                                = arg_int1(NULL,"new_k", NULL,"number of blocks in repartition");
  struct arg_lit *diffusive_refinement                 = arg_lit0(NULL, "diffusive_refinement", "choice of diffusive refinement");
  struct arg_lit *trunccons_limited                    = arg_lit0(NULL, "trunccons_limited", "choice of limited trunccons");
  struct arg_lit *trunccons_limited_useAllNeighbors    = arg_lit0(NULL, "trunccons_limited_useAllNeighbors", "choice of using all neighbors");
  struct arg_lit *scale_balancing                      = arg_lit0(NULL, "scale_balancing", "choice of scale balancing");
  struct arg_lit *flow_based_balancing                 = arg_lit0(NULL, "flow_based_balancing", "choice of flow based balancing");
  struct arg_lit *greedy_balancing                     = arg_lit0(NULL, "greedy_balancing", "choice of greedy balancing");
  struct arg_int *trunccons_numDiffIters               = arg_int0(NULL,"trunccons_numDiffIters", NULL, "number of diffusion iterations of trunccons");
  struct arg_int *trunccons_numConsols                 = arg_int0(NULL,"trunccons_numConsols", NULL,"number of consolidations in trunccons");
  struct arg_int *trunccons_limited_sizeOfNeighborhood = arg_int0(NULL,"trunccons_limited_sizeOfNeighborhood", NULL,"size of neighborhood in trunccons_limited");
  struct arg_int *scale_balancing_iters                = arg_int0(NULL,"scale_balancing_iters", NULL,"number of iterations in scale-balancing");
#endif
  
#ifdef MODE_MAP
  struct arg_str *mapping_algo                    = arg_str0(NULL, "mapping_algo", NULL, "algorithm used for mapping, default: initial");
  struct arg_int *numberHierarchies               = arg_int0(NULL, "numberHierarchies", NULL, "Number of hierachies used by Alma.");
  struct arg_str *proc_graph_type                 = arg_str0(NULL, "proc_graph_type", NULL, "type of processor graph");
  struct arg_str *proc_graph_extensions           = arg_str0(NULL, "proc_graph_extensions", NULL, "e.g., '16x16' for 2D grid of '4x3x3' for tree, depth 3, top-down");
  struct arg_str *proc_graph_bandwidths           = arg_str0(NULL, "proc_graph_bandwidths", NULL, "e.g., 3x5x7x9 for 16x16 torus or tree of depth 4, top-down");
#endif
  
  struct arg_end *end                                  = arg_end(100);
  
  // Define argtable.
  void* argtable[] = {
    help, filename, user_seed,
#ifdef MODE_DEVEL
    k, graph_weighted, imbalance, edge_rating_tiebreaking, 
    matching_type, edge_rating, rate_first_level_inner_outer, first_level_random_matching, 
    aggressive_random_levels, gpa_grow_internal, match_islands, stop_rule, num_vert_stop_factor,
    initial_partition, initial_partitioning_repetitions, disable_refined_bubbling, 
    bubbling_iterations, initial_partition_optimize, bipartition_post_fm_limit, bipartition_post_ml_limit, bipartition_tries,
    bipartition_algorithm,
    permutation_quality, permutation_during_refinement, enforce_balance,
    refinement_scheduling_algorithm, bank_account_factor, refinement_type, 
    fm_search_limit, flow_region_factor, most_balanced_flows,toposort_iterations, 
    kway_rounds, kway_search_stop_rule, kway_fm_limits, kway_adaptive_limits_alpha, 
    enable_corner_refinement, disable_qgraph_refinement,local_multitry_fm_alpha, local_multitry_rounds,
    global_cycle_iterations, use_wcycles, wcycle_no_new_initial_partitioning, use_fullmultigrid, use_vcycle,level_split, 
    enable_convergence, compute_vertex_separator, suppress_output, 
    input_partition, preconfiguration, only_first_level, disable_max_vertex_weight_constraint, 
    recursive_bipartitioning, use_bucket_queues, time_limit, unsuccessful_reps, local_partitioning_repetitions, 
    mh_pool_size, mh_plain_repetitions, mh_disable_nc_combine, mh_disable_cross_combine, mh_enable_tournament_selection,       
    mh_disable_combine, mh_enable_quickstart, mh_disable_diversify_islands, mh_flip_coin, mh_initial_population_fraction, 
    mh_num_ncs_to_compute, mh_print_log,mh_sequential_mode, mh_optimize_communication_volume, mh_enable_tabu_search,
    mh_disable_diversify, mh_diversify_best, mh_cross_combine_original_k, disable_balance_singletons, initial_partition_optimize_fm_limits,
    initial_partition_optimize_multitry_fm_alpha, initial_partition_optimize_multitry_rounds,
    enable_omp, 
    amg_iterations,
    kaba_neg_cycle_algorithm, kabaE_internal_bal, kaba_internal_no_aug_steps_aug, 
    kaba_packing_iterations, kaba_flip_packings, kaba_lsearch_p, kaffpa_perfectly_balanced_refinement, 
    kaba_unsucc_iterations, kaba_disable_zero_weight_cycles,
    maxT, maxIter, minipreps, mh_penalty_for_unconnected, mh_enable_kabapE,

#elif defined MODE_KAFFPA
    k, imbalance,  
    preconfiguration, 
    time_limit, 
    input_partition,
    enforce_balance, 
    filename_output, 

#elif defined MODE_PARTITION_INFO
    k,
    input_partition,

#elif defined MODE_PART
    k, imbalance,  
    preconfiguration,
    edge_rating, 
    time_limit, 
    input_partition,
    enforce_balance, 
    filename_output, 

#elif defined MODE_REPART
    k, imbalance,  
    preconfiguration,
    edge_rating,
    time_limit, 
    preconfiguration,
    input_partition,
    enforce_balance, 
    filename_output,
    filename_old_graph,
    filename_ancestors,
    new_imbalance,
    new_k,
    diffusive_refinement,
    trunccons_limited, trunccons_limited_useAllNeighbors,
    scale_balancing, flow_based_balancing, greedy_balancing,
    trunccons_numDiffIters, trunccons_numConsols, trunccons_limited_sizeOfNeighborhood,
    scale_balancing_iters,
    
#elif defined MODE_MAP
    k, imbalance,
    preconfiguration,
    edge_rating,
    time_limit, 
    preconfiguration,
    input_partition,
    enforce_balance, 
    filename_output,
    mapping_algo,
    numberHierarchies,
    multikahip_preconfiguration,
    proc_graph_type,
    proc_graph_extensions,
    proc_graph_bandwidths,

#elif defined MODE_PARTITIONTOVERTEXSEPARATOR
    k, input_partition, 
    filename_output, 
    
#elif defined MODE_KAFFPAE
    k, imbalance, 
    preconfiguration,  
    time_limit,  
    mh_enable_quickstart, 
    mh_print_log, mh_optimize_communication_volume, 
    mh_enable_tabu_search,
    maxT, maxIter,  
    mh_enable_kabapE,
    kabaE_internal_bal,  
    input_partition,
    filename_output, 
#elif defined MODE_LABELPROPAGATION
    cluster_upperbound,
    label_propagation_iterations,
    filename_output, 
#endif
    end
  };
  // Parse arguments.
  int nerrors = arg_parse(argn, argv, argtable);
  
  // Catch case that help was requested.
  if (help->count > 0) {
    printf("Usage: %s", progname);
    arg_print_syntax(stdout, argtable, "\n");
    arg_print_glossary(stdout, argtable,"  %-40s %s\n");
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return 1;
  }

  
  if (nerrors > 0) {
    arg_print_errors(stderr, end, progname);
    printf("Try '%s --help' for more information.\n",progname);
    //arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return 1; 
  }
  
  if (k->count > 0) {
    partition_config.k = k->ival[0];
  }

  if(filename->count > 0) {
    graph_filename = filename->sval[0];
  }
  
  recursive = false;
  
  configuration cfg;
    cfg.standard(partition_config);
    
#ifdef MODE_KAFFPA
    cfg.eco(partition_config);
#endif
    
#ifdef MODE_PART
    cfg.eco(partition_config);
#endif
    
#ifdef MODE_REPART
    cfg.eco(partition_config);
#endif
    
#ifdef MODE_MAP
    cfg.eco(partition_config);
#endif
    
#ifdef MODE_KAFFPAE
    cfg.strong(partition_config);
#endif
    
#ifdef MODE_PARTITIONTOVERTEXSEPARATOR
    cfg.strong(partition_config);
#endif
    
#ifdef MODE_LABELPROPAGATION
    cfg.strong(partition_config);
#endif
    
    if(preconfiguration->count > 0) {
      if(strcmp("strong", preconfiguration->sval[0]) == 0) {
	cfg.strong(partition_config);
    } else if (strcmp("eco", preconfiguration->sval[0]) == 0) {
	cfg.eco(partition_config);
      } else if (strcmp("fast", preconfiguration->sval[0]) == 0) {
	cfg.fast(partition_config);
      } else if (strcmp("fastsocial", preconfiguration->sval[0]) == 0) {
	cfg.fastsocial(partition_config);
      } else if (strcmp("ecosocial", preconfiguration->sval[0]) == 0) {
	cfg.ecosocial(partition_config);
      } else if (strcmp("strongsocial", preconfiguration->sval[0]) == 0) {
	cfg.strongsocial(partition_config);
      } else {
	fprintf(stderr, "Invalid preconfiguration variant: \"%s\"\n", preconfiguration->sval[0]);
      exit(0);
      }
    }


    if(filename_output->count > 0) {
      partition_config.filename_output = filename_output->sval[0];
    } else {
      partition_config.filename_output = "";
    }

    if(initial_partition_optimize->count > 0) {
      partition_config.initial_partition_optimize = true;
    }
    
    if(disable_balance_singletons->count > 0) {
      partition_config.use_balance_singletons = false;
    }
    
    if(mh_disable_nc_combine->count > 0) {
      partition_config.mh_disable_nc_combine = true;
    }
    
    if(mh_disable_cross_combine->count > 0) {
      partition_config.mh_disable_cross_combine = true;
    }
    
    if( imbalance->count > 0) {
      partition_config.epsilon = imbalance->dval[0];
    }
    
    if(mh_disable_combine->count > 0) {
      partition_config.mh_disable_combine = true;
    }
    
    if(mh_optimize_communication_volume->count > 0) {
      partition_config.mh_optimize_communication_volume = true;
    }
    
    if(mh_enable_tournament_selection->count > 0) {
      partition_config.mh_enable_tournament_selection = true;
    }
    
    if(amg_iterations->count > 0) {
      partition_config.amg_iterations = amg_iterations->ival[0];
    }
    
    if(kabaE_internal_bal->count > 0) {
      partition_config.kabaE_internal_bal = kabaE_internal_bal->dval[0];
    }
    
    if(kaba_internal_no_aug_steps_aug->count > 0) {
      partition_config.kaba_internal_no_aug_steps_aug = kaba_internal_no_aug_steps_aug->ival[0];
    }
    
    if(kaffpa_perfectly_balanced_refinement->count > 0) {
      partition_config.kaffpa_perfectly_balanced_refinement = true;
    }

    if(kaba_disable_zero_weight_cycles->count > 0) {
      partition_config.kaba_enable_zero_weight_cycles = false;
    }
    
    if(kaba_unsucc_iterations->count > 0) {
      partition_config.kaba_unsucc_iterations = kaba_unsucc_iterations->ival[0];
    }

    if(kaba_flip_packings->count > 0) {
      partition_config.kaba_flip_packings = true;
    }
    
    if (kaba_lsearch_p->count) {
      if(strcmp("coindiff", kaba_lsearch_p->sval[0]) == 0) {
	partition_config.kaba_lsearch_p = COIN_DIFFTIE;
      } else if (strcmp("nocoindiff",kaba_lsearch_p->sval[0]) == 0) {
	partition_config.kaba_lsearch_p = NOCOIN_DIFFTIE;
      } else if (strcmp("coinrnd", kaba_lsearch_p->sval[0]) == 0) {
	partition_config.kaba_lsearch_p = COIN_RNDTIE;
      } else if (strcmp("nocoinrnd", kaba_lsearch_p->sval[0]) == 0) {
	partition_config.kaba_lsearch_p = NOCOIN_RNDTIE;
      } else {
	fprintf(stderr, "Invalid combine variant: \"%s\"\n", kaba_lsearch_p->sval[0]);
	exit(0);
      }
    }
    
    if(maxT->count > 0) {
      partition_config.maxT = maxT->ival[0];
    }
    
    if(maxIter->count > 0) {
      partition_config.maxIter = maxIter->ival[0];
    }
    
    
    if(mh_enable_tabu_search->count > 0) {
      partition_config.mh_enable_gal_combine = true;
    }
    
    if(kaba_packing_iterations->count > 0) {
      partition_config.kaba_packing_iterations = kaba_packing_iterations->ival[0];
    }
    
    if(mh_flip_coin->count > 0) {
      partition_config.mh_flip_coin = mh_flip_coin->ival[0];
    }
    
    if(mh_initial_population_fraction->count > 0) {
      partition_config.mh_initial_population_fraction = mh_initial_population_fraction->ival[0];
    }
    
    if(minipreps->count > 0) {
      partition_config.minipreps = minipreps->ival[0];
    }
    
    
    if(mh_enable_quickstart->count > 0) {
      partition_config.mh_enable_quickstart = true;
    }
    
    if(mh_disable_diversify_islands->count > 0) {
      partition_config.mh_disable_diversify_islands = true;
    }
    
    if(gpa_grow_internal->count > 0) {
      partition_config.gpa_grow_paths_between_blocks = false;
    }

    if(suppress_output->count > 0) {
      suppress_program_output = true;
    }

    if(mh_print_log->count > 0) {
      partition_config.mh_print_log = true;
    }

    if(use_bucket_queues->count > 0) {
      partition_config.use_bucket_queues = true;
    }

    if(recursive_bipartitioning->count > 0 ) {
      recursive = true;
    }

    if(time_limit->count > 0) {
      partition_config.time_limit = time_limit->dval[0];
    }

    if(unsuccessful_reps->count > 0) {
      partition_config.no_unsuc_reps = unsuccessful_reps->ival[0];
    }

    if(mh_pool_size->count > 0) {
      partition_config.mh_pool_size = mh_pool_size->ival[0];
    }

    if(mh_penalty_for_unconnected->count > 0) {
      partition_config.mh_penalty_for_unconnected = true;
    }

    if(mh_enable_kabapE->count > 0) {
      partition_config.kabapE = true;
    }

    if(initial_partition_optimize_multitry_fm_alpha->count > 0) {
      partition_config.initial_partition_optimize_multitry_fm_alpha = initial_partition_optimize_multitry_fm_alpha->ival[0];
    }

    if(initial_partition_optimize_multitry_rounds->count > 0) {
      partition_config.initial_partition_optimize_multitry_rounds = initial_partition_optimize_multitry_rounds->ival[0];
    }

    if(initial_partition_optimize_fm_limits->count > 0) {
      partition_config.initial_partition_optimize_fm_limits = initial_partition_optimize_fm_limits->ival[0];
    }

    if(mh_disable_diversify->count > 0) {
      partition_config.mh_diversify = false;
    }

    if(mh_diversify_best->count > 0) {
      partition_config.mh_diversify_best = true;
    }

    if(enforce_balance->count > 0) {
      partition_config.kaffpa_perfectly_balance = true;
    }

    if(mh_plain_repetitions->count > 0) {
      partition_config.mh_plain_repetitions = true;
    }

    if(local_partitioning_repetitions->count > 0) {
      partition_config.local_partitioning_repetitions = local_partitioning_repetitions->ival[0];
    }

    if(only_first_level->count > 0) {
      partition_config.only_first_level = true;
    }

    if(mh_cross_combine_original_k->count > 0) {
      partition_config.mh_cross_combine_original_k = true;
    }

    if(mh_sequential_mode->count > 0) {
      partition_config.mh_no_mh = true;
    }

    if(enable_omp->count > 0) {
      partition_config.enable_omp = true;
    }

    if(compute_vertex_separator->count > 0) {
      partition_config.compute_vertex_separator = true;
    }

    if(most_balanced_flows->count > 0) {
      partition_config.most_balanced_minimum_cuts = true;
    }

    if(use_wcycles->count > 0) {
      partition_config.use_wcycles = true; 
    }

    if(enable_convergence->count > 0) {
      partition_config.no_change_convergence = true; 
    }

    if(use_fullmultigrid->count > 0) {
      partition_config.use_fullmultigrid = true; 
    }

    if(use_vcycle->count > 0) {
      partition_config.use_fullmultigrid = false; 
      partition_config.use_wcycles       = false; 
    }


    if(toposort_iterations->count > 0) {
      partition_config.toposort_iterations = toposort_iterations->ival[0]; 
    }

    if(bipartition_tries->count > 0) {
      partition_config.bipartition_tries = bipartition_tries->ival[0]; 
    }

    if(bipartition_post_fm_limit->count > 0) {
      partition_config.bipartition_post_fm_limits = bipartition_post_fm_limit->ival[0]; 
    }

    if(bipartition_post_ml_limit->count > 0) {
      partition_config.bipartition_post_ml_limits = bipartition_post_ml_limit->ival[0]; 
    }

    if(disable_max_vertex_weight_constraint->count > 0) {
      partition_config.disable_max_vertex_weight_constraint = true;
    }

    if(num_vert_stop_factor->count > 0) {
      partition_config.num_vert_stop_factor = num_vert_stop_factor->ival[0]; 
    }

    if(local_multitry_rounds->count > 0) {
      partition_config.local_multitry_rounds = local_multitry_rounds->ival[0]; 
    }

    if(local_multitry_fm_alpha->count > 0) {
      partition_config.local_multitry_fm_alpha = local_multitry_fm_alpha->ival[0]; 
    }

    if(wcycle_no_new_initial_partitioning->count > 0) {
      partition_config.no_new_initial_partitioning = true; 
    }

    if(graph_weighted->count > 0) {
      is_graph_weighted = true;                
    }

    if(disable_refined_bubbling->count > 0) {
      partition_config.refined_bubbling = false;
    }

    if(input_partition->count > 0) {
      partition_config.input_partition = input_partition->sval[0];
    } else {
      partition_config.input_partition = "";
    }

    if(global_cycle_iterations->count > 0) {
      partition_config.global_cycle_iterations = global_cycle_iterations->ival[0];
    }

    if(level_split->count > 0) {
      partition_config.level_split = level_split->ival[0];
    }

    if(disable_qgraph_refinement->count > 0) {
      partition_config.quotient_graph_refinement_disabled = true;
    }

    if(bubbling_iterations->count>0) {
      partition_config.bubbling_iterations = bubbling_iterations->ival[0];
    }

    if(kway_fm_limits->count>0) {
      partition_config.kway_fm_search_limit = kway_fm_limits->ival[0];
    }

    if(kway_rounds->count>0) {
      partition_config.kway_rounds = kway_rounds->ival[0];
    }

    if(enable_corner_refinement->count > 0) {
      partition_config.corner_refinement_enabled = true;
    }

    if(match_islands->count > 0) {
      partition_config.match_islands = true;
    }

    if(aggressive_random_levels->count > 0) {
      partition_config.aggressive_random_levels = aggressive_random_levels->ival[0];
    }

    if(rate_first_level_inner_outer->count > 0) {
      partition_config.rate_first_level_inner_outer = true;
    }

    if(user_seed->count > 0) {
      partition_config.seed = user_seed->ival[0];
    }

    if(fm_search_limit->count > 0) {
      partition_config.fm_search_limit = fm_search_limit->ival[0];
    }

    if(bank_account_factor->count > 0) {
      partition_config.bank_account_factor = bank_account_factor->dval[0];
    }

    if(flow_region_factor->count > 0) {
      partition_config.flow_region_factor = flow_region_factor->dval[0];
    }

    if(kway_adaptive_limits_alpha->count > 0) {
      partition_config.kway_adaptive_limits_alpha = kway_adaptive_limits_alpha->dval[0];
    }

    if(imbalance->count > 0) {
      partition_config.imbalance = imbalance->dval[0];
    }

    if(initial_partitioning_repetitions->count > 0) {
      partition_config.initial_partitioning_repetitions = initial_partitioning_repetitions->ival[0];
    }

    if(edge_rating_tiebreaking->count > 0) {
      partition_config.edge_rating_tiebreaking = true;
    }

    if(first_level_random_matching->count > 0) {
      partition_config.first_level_random_matching = true;
    } else {
      partition_config.first_level_random_matching = false;
    }

    if(kaba_neg_cycle_algorithm->count > 0) {
      if (strcmp("playfield",kaba_neg_cycle_algorithm->sval[0]) == 0) {
	partition_config.cycle_refinement_algorithm = CYCLE_REFINEMENT_ALGORITHM_PLAYFIELD;
      } else if (strcmp("ultramodel",kaba_neg_cycle_algorithm->sval[0]) == 0) {
	partition_config.cycle_refinement_algorithm = CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL;
      } else if (strcmp("ultramodelplus",kaba_neg_cycle_algorithm->sval[0]) == 0) {
	partition_config.cycle_refinement_algorithm = CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL_PLUS;
      } else {
	fprintf(stderr, "Invalid balanced refinement operator: \"%s\"\n", kaba_neg_cycle_algorithm->sval[0]);
	exit(0);
      }
    }

  
    if (edge_rating->count > 0) {
      if(strcmp("expansionstar", edge_rating->sval[0]) == 0) {
	partition_config.edge_rating = EXPANSIONSTAR;
      } else if (strcmp("expansionstar2", edge_rating->sval[0]) == 0) {
	partition_config.edge_rating = EXPANSIONSTAR2;
      } else if (strcmp("expansionstar2algdist", edge_rating->sval[0]) == 0) {
	partition_config.edge_rating = EXPANSIONSTAR2ALGDIST;
      } else if (strcmp("geom", edge_rating->sval[0]) == 0) {
	partition_config.edge_rating = PSEUDOGEOM;
      } else if (strcmp("resistance", edge_rating->sval[0]) == 0) {
	partition_config.edge_rating = RESISTANCE;
      } else if (strcmp("zero", edge_rating->sval[0]) == 0) {
	partition_config.edge_rating = ZERO;
      } else {
	fprintf(stderr, "Invalid edge rating variant: \"%s\"\n", edge_rating->sval[0]);
	exit(0);
      }
    }
    
    if(bipartition_algorithm->count > 0) {
      if(strcmp("bfs", bipartition_algorithm->sval[0]) == 0) {
	partition_config.bipartition_algorithm = BIPARTITION_BFS;
      } else if (strcmp("fm", bipartition_algorithm->sval[0]) == 0) {
	partition_config.bipartition_algorithm = BIPARTITION_FM;
      } else {
	fprintf(stderr, "Invalid bipartition algorthim: \"%s\"\n", bipartition_algorithm->sval[0]);
	exit(0);
      }
    }

    if(refinement_scheduling_algorithm->count > 0) {
      if(strcmp("fast", refinement_scheduling_algorithm->sval[0]) == 0) {
	partition_config.refinement_scheduling_algorithm = REFINEMENT_SCHEDULING_FAST;
      } else if (strcmp("active_blocks", refinement_scheduling_algorithm->sval[0]) == 0) {
	partition_config.refinement_scheduling_algorithm = REFINEMENT_SCHEDULING_ACTIVE_BLOCKS;
      } else if (strcmp("active_blocks_kway", refinement_scheduling_algorithm->sval[0]) == 0) {
	partition_config.refinement_scheduling_algorithm = REFINEMENT_SCHEDULING_ACTIVE_BLOCKS_REF_KWAY;
      } else {
	fprintf(stderr, "Invalid refinement scheduling variant: \"%s\"\n", refinement_scheduling_algorithm->sval[0]);
	exit(0);
      }
    }

    if(stop_rule->count > 0) {
      if(strcmp("simple", stop_rule->sval[0]) == 0) {
	partition_config.stop_rule = STOP_RULE_SIMPLE;
      } else if (strcmp("multiplek", stop_rule->sval[0]) == 0) {
	partition_config.stop_rule = STOP_RULE_MULTIPLE_K;
      } else if (strcmp("strong", stop_rule->sval[0]) == 0) {
	partition_config.stop_rule = STOP_RULE_STRONG;
      } else {
	fprintf(stderr, "Invalid stop rule: \"%s\"\n", stop_rule->sval[0]);
	exit(0);
      }
    }

    if(kway_search_stop_rule->count > 0) {
      if(strcmp("simple", kway_search_stop_rule->sval[0]) == 0) {
	partition_config.kway_stop_rule = KWAY_SIMPLE_STOP_RULE;
      } else if (strcmp("adaptive", kway_search_stop_rule->sval[0]) == 0) {
	partition_config.kway_stop_rule = KWAY_ADAPTIVE_STOP_RULE;
      } else {
	fprintf(stderr, "Invalid kway stop rule: \"%s\"\n", kway_search_stop_rule->sval[0]);
	exit(0);
      }
    }

    if(permutation_quality->count > 0) {
      if(strcmp("none", permutation_quality->sval[0]) == 0) {
	partition_config.permutation_quality = PERMUTATION_QUALITY_NONE;
      } else if (strcmp("fast", permutation_quality->sval[0]) == 0) {
	partition_config.permutation_quality = PERMUTATION_QUALITY_FAST;
      } else if (strcmp("good", permutation_quality->sval[0]) == 0) {
	partition_config.permutation_quality = PERMUTATION_QUALITY_GOOD;
      } else {
	fprintf(stderr, "Invalid permutation quality variant: \"%s\"\n", permutation_quality->sval[0]);
	exit(0);
      }

    }

    if(permutation_during_refinement->count > 0) {
      if(strcmp("none", permutation_during_refinement->sval[0]) == 0) {
	partition_config.permutation_during_refinement = PERMUTATION_QUALITY_NONE;
      } else if (strcmp("fast", permutation_during_refinement->sval[0]) == 0) {
	partition_config.permutation_during_refinement = PERMUTATION_QUALITY_FAST;
      } else if (strcmp("good", permutation_during_refinement->sval[0]) == 0) {
	partition_config.permutation_during_refinement = PERMUTATION_QUALITY_GOOD;
      } else {
	fprintf(stderr, "Invalid permutation quality during refinement variant: \"%s\"\n", permutation_during_refinement->sval[0]);
	exit(0);
      }
    }

    if(matching_type->count > 0) {
      if(strcmp("random", matching_type->sval[0]) == 0) {
	partition_config.matching_type = MATCHING_RANDOM;
      } else if (strcmp("gpa", matching_type->sval[0]) == 0) {
	partition_config.matching_type = MATCHING_GPA;
      } else if (strcmp("randomgpa", matching_type->sval[0]) == 0) {
	partition_config.matching_type = MATCHING_RANDOM_GPA;
      } else {
	fprintf(stderr, "Invalid matching variant: \"%s\"\n", matching_type->sval[0]);
	exit(0);
      }
    }

    if(refinement_type->count > 0) {
      if(strcmp("fm", refinement_type->sval[0]) == 0) {
	partition_config.refinement_type = REFINEMENT_TYPE_FM;
      } else if (strcmp("fm_flow", refinement_type->sval[0]) == 0) {
	partition_config.refinement_type = REFINEMENT_TYPE_FM_FLOW;
      } else if (strcmp("flow", refinement_type->sval[0]) == 0) {
	partition_config.refinement_type = REFINEMENT_TYPE_FLOW;
      } else {
	fprintf(stderr, "Invalid refinement type variant: \"%s\"\n", refinement_type->sval[0]);
	exit(0);
      }
    }

    if(initial_partition->count > 0) { 
      if (strcmp("recursive", initial_partition->sval[0]) == 0) {
	partition_config.initial_partitioning_type = INITIAL_PARTITIONING_RECPARTITION;
      } else {
	fprintf(stderr, "Invalid initial partition variant: \"%s\"\n", initial_partition->sval[0]);
	exit(0);
      }
    }

    if(label_propagation_iterations->count > 0) {
      partition_config.label_iterations = label_propagation_iterations->ival[0];
    }
	
    if(cluster_upperbound->count > 0) {
      partition_config.cluster_upperbound = cluster_upperbound->ival[0];
    } else {
      partition_config.cluster_upperbound = std::numeric_limits< NodeWeight >::max()/2;
    }

#ifdef MODE_REPART
    if(filename_old_graph->count > 0) {
      partition_config.filename_old_graph = filename_old_graph->sval[0];
    } else {
      partition_config.filename_old_graph = "";
    }
  
    if(filename_ancestors->count > 0) {
      partition_config.filename_ancestors = filename_ancestors->sval[0];
    } else {
      partition_config.filename_ancestors = "";
    }
  
    if(new_imbalance->count > 0) {
      partition_config.new_imbalance = new_imbalance->dval[0];
    } else {
      partition_config.new_imbalance = 3;
    }
  
    if(new_k->count > 0) {
      partition_config.new_k = new_k->ival[0];
    } else {
      partition_config.new_k = 0.03;
    }
  
    if(diffusive_refinement->count > 0) {
      partition_config.diffusive_refinement = true;
    }
  
    if(trunccons_limited->count > 0) {
      partition_config.trunccons_limited = true;
    }
  
    if(trunccons_limited_useAllNeighbors->count > 0) {
      partition_config.trunccons_limited_useAllNeighbors = true;
    }
  
    if(scale_balancing->count > 0) {
      partition_config.scale_balancing = true;
    }
  
    if(flow_based_balancing->count > 0) {
      partition_config.flow_based_balancing = true;
    }

    if(greedy_balancing->count > 0) {
      partition_config.greedy_balancing = true;
    }
  
    if(trunccons_numDiffIters->count > 0) {
      partition_config.trunccons_numDiffIters = trunccons_numDiffIters->ival[0];
    } else {
      partition_config.trunccons_numDiffIters = 14;
    }
  
    if(trunccons_numConsols->count > 0) {
      partition_config.trunccons_numConsols = trunccons_numConsols->ival[0];
    } else {
      partition_config.trunccons_numConsols = 9;
    }
  
    if(trunccons_limited_sizeOfNeighborhood->count > 0) {
      partition_config.trunccons_limited_sizeOfNeighborhood = trunccons_limited_sizeOfNeighborhood->ival[0];
    } else {
      partition_config.trunccons_limited_sizeOfNeighborhood = 8;
    }
  
    if(scale_balancing_iters->count > 0) {
      partition_config.scale_balancing_iters = scale_balancing_iters->ival[0];
    } else {
      partition_config.scale_balancing_iters = 16;
    }
#endif
  
#ifdef MODE_MAP
    if(mapping_algo->count > 0) { 
      if(strcmp("initial", mapping_algo->sval[0]) == 0) {
	partition_config.mapping_algo = mapping_algo->sval[0];
      } else if (strcmp("random", mapping_algo->sval[0]) == 0) {
	partition_config.mapping_algo = mapping_algo->sval[0];
      } else if (strcmp("drb", mapping_algo->sval[0]) == 0) {
	partition_config.mapping_algo = mapping_algo->sval[0];
      } else if (strcmp("rcm", mapping_algo->sval[0]) == 0) {
	partition_config.mapping_algo = mapping_algo->sval[0];
      } else if (strcmp("multisimple", mapping_algo->sval[0]) == 0) {
	partition_config.mapping_algo = mapping_algo->sval[0];
      } else if (strcmp("multikahip", mapping_algo->sval[0]) == 0) {
	partition_config.mapping_algo = mapping_algo->sval[0];
      } else if (strcmp("brandfass", mapping_algo->sval[0]) == 0) {
	partition_config.mapping_algo = mapping_algo->sval[0];
      } else if (strcmp("brandfassL", mapping_algo->sval[0]) == 0) {
	partition_config.mapping_algo = mapping_algo->sval[0];
      } else if (strcmp("hoefler", mapping_algo->sval[0]) == 0) {
	partition_config.mapping_algo = mapping_algo->sval[0];
      } else if (strcmp("hoeflerL", mapping_algo->sval[0]) == 0) {
	partition_config.mapping_algo = mapping_algo->sval[0];
      } else if (strcmp("alma", mapping_algo->sval[0]) == 0) { /* maria */
      	partition_config.mapping_algo = mapping_algo->sval[0];
      } else {
	fprintf(stderr, "Invalid mapping algorithm: \"%s\"\n", mapping_algo->sval[0]);
	exit(0);
      }
    }

    if(numberHierarchies->count > 0) { 
      partition_config.numberHierarchies = numberHierarchies->ival[0];
    } else {
      partition_config.numberHierarchies = 0;
    }

    if(proc_graph_type->count > 0) {
      if(strcmp("grid", proc_graph_type->sval[0]) == 0) {
	partition_config.proc_graph_type = proc_graph_type->sval[0];
      } else if (strcmp("torus", proc_graph_type->sval[0]) == 0) {
	partition_config.proc_graph_type = proc_graph_type->sval[0];
      } else if (strcmp("tree", proc_graph_type->sval[0]) == 0) {
	partition_config.proc_graph_type = proc_graph_type->sval[0];
      } else if (strcmp("hypercube", proc_graph_type->sval[0]) == 0) {
	partition_config.proc_graph_type = proc_graph_type->sval[0];
      } else {
	fprintf(stderr, "Invalid processor graph type: \"%s\"\n", proc_graph_type->sval[0]);
	exit(0);
      }
    }
  
    if(proc_graph_extensions->count > 0) {
      partition_config.proc_graph_extensions = proc_graph_extensions->sval[0];
    } else {
      partition_config.proc_graph_extensions = "";
    }
  
    if(proc_graph_bandwidths->count > 0) {
      partition_config.proc_graph_bandwidths = proc_graph_bandwidths->sval[0];
    } else {
      partition_config.proc_graph_bandwidths = "";
    }
  
    if(multikahip_preconfiguration->count > 0) {
      if(strcmp("strong", multikahip_preconfiguration->sval[0]) == 0) {
	partition_config.multikahip_preconfiguration = multikahip_preconfiguration->sval[0];
      } else if (strcmp("eco", multikahip_preconfiguration->sval[0]) == 0) {
	partition_config.multikahip_preconfiguration = multikahip_preconfiguration->sval[0];
      } else if (strcmp("fast", multikahip_preconfiguration->sval[0]) == 0) {
	partition_config.multikahip_preconfiguration = multikahip_preconfiguration->sval[0];
      } else if (strcmp("fastsocial", multikahip_preconfiguration->sval[0]) == 0) {
	partition_config.multikahip_preconfiguration = multikahip_preconfiguration->sval[0];
      } else if (strcmp("ecosocial", multikahip_preconfiguration->sval[0]) == 0) {
	partition_config.multikahip_preconfiguration = multikahip_preconfiguration->sval[0];
      } else if (strcmp("strongsocial", multikahip_preconfiguration->sval[0]) == 0) {
	partition_config.multikahip_preconfiguration = multikahip_preconfiguration->sval[0];
      } else {
	fprintf(stderr, "Invalid multikahip_preconfiguration variant: \"%s\"\n", multikahip_preconfiguration->sval[0]);
	exit(0);
      }
    }

#endif

    return 0;
}

/***********************************************************************************************************/
/******************************************** parse_repart_parameters() ************************************/
/***********************************************************************************************************/
int parse_repart_parameters(int argn, char **argv,
			    PartitionConfig & partition_config,
			    PartitionConfig & repartition_config,
			    std::string & graph_filename, 
			    bool & is_graph_weighted, 
			    bool & suppress_program_output, 
			    bool & recursive) {
  int ret=parse_parameters(argn,argv, partition_config, graph_filename, is_graph_weighted, suppress_program_output, recursive);
  if(ret == 0) {
    repartition_config = partition_config;
    repartition_config.imbalance = partition_config.new_imbalance;
    repartition_config.k = partition_config.new_k;
    repartition_config.input_partition = "";
  }
  return ret;
}
  
/***********************************************************************************************************/
/********************************************* parse_map_parameters() **************************************/
/***********************************************************************************************************/
int parse_map_parameters(int argn, char **argv,
			 PartitionConfig & partition_config,
			 PartitionConfig & map_config,
			 std::string & graph_filename,
			 bool & is_graph_weighted,
			 bool & suppress_program_output,
			 bool & recursive) {
  int ret = parse_parameters(argn,argv, partition_config, graph_filename, is_graph_weighted, suppress_program_output, recursive);
  printf("parse_parameters:::ret is %d\n", ret);
  if(ret == 0) {
    map_config = partition_config;    
    configuration cfg;
    if(map_config.multikahip_preconfiguration == "fast") {
      cfg.fast(map_config);
    } else {
      if((map_config.multikahip_preconfiguration == "eco") || (map_config.multikahip_preconfiguration == "")) {
	cfg.eco(map_config);
      } else {
	if(map_config.multikahip_preconfiguration == "strong") {
	  cfg.strong(map_config);
	} else {
	  if(map_config.multikahip_preconfiguration == "fastsocial") {
	    cfg.fastsocial(map_config);
	  } else {
	    if(map_config.multikahip_preconfiguration == "ecosocial") {
	      cfg.ecosocial(map_config);
	    } else {
	      if(map_config.multikahip_preconfiguration == "strongsocial") {
		cfg.strongsocial(map_config);
	      }
	    }
	  }
	}
      }
    }
  }
  return ret;
}

#endif /* end of include guard: PARSE_PARAMETERS_GPJMGSM8 */
