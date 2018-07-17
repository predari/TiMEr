#include "data_structure/graph_access.h"
#include "../lib/data_structure/graph_access.h"
#include "./errors.h"
#include "./processor_graph_access.h"
#include <vector>
#include <algorithm>

using namespace std;

double random01(void);
int getPower(int size);
void build2DTreeRec(int start, int total, int size, vector<int>* tree);
void build3DTreeRec(int start, int total, int size, vector<int>* tree);
vector< vector<int> > buildLevDist(int level);
int elementOfWeightVector(int size, int dest, int target);
void allPathsShortestPaths (vector< vector<int> >* predecessorMatrix, vector< vector<int> >* distanceMatrix, graph_access* graph);
void build2DGrid(processor_graph_access & proc, int size, vector<int> & bandwidths);
void build3DGrid (processor_graph_access & proc, int size, vector<int> & bandwidths);
void build2DTorus (processor_graph_access & proc, int size, vector<int> & bandwidths);
void build3DTorus (processor_graph_access & proc, int size, vector<int> & bandwidths);
void buildFatTree(processor_graph_access & proc, vector<int> & dimensions, vector<int> & bandwidths);
void buildHypercube(processor_graph_access & proc, int dimension);
