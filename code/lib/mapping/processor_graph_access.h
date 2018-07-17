/* Roland Glantz <roland.glantz@kit.edu>
 * Adapted from GraphCreation.h by Alexander Noe
 */

#ifndef PROCESSOR_GRAPH_H_
#define PROCESSOR_GRAPH_H_

using namespace std;
#include "../lib/data_structure/graph_access.h"

#include <vector>

class processor_graph_access : public graph_access
{
 protected:

  vector<int64_t> hl;//Hamming labels on vertices of processor graph

 public:
  bool refine;
  int dim;
  vector< vector<int> > distanceMatrix;
  vector< vector<int> > predecessorMatrix;
  vector<int> splitDeg;
  vector<int> treeSort;
  vector< vector<int> > levelDist;

  processor_graph_access();

  ~processor_graph_access();

  inline void resize_HammingLabels(NodeID size) {
    hl.resize(size);
  }

  inline void setHammingLabel(NodeID v, int64_t l) {
    hl[v] = l;
  }

  inline int64_t getHammingLabel(NodeID v) {
    return(hl[v]);
  }

  void getDist2v(processor_graph_access & P, vector<int64_t> &dist2v, int64_t v);

  bool setHammingLabels(processor_graph_access & P);

};

processor_graph_access createFatTree (char* size, char* weights);
processor_graph_access create2DTorus (int size, char* weights);
processor_graph_access create3DTorus (int size, char* weights);
processor_graph_access create2DGrid (int size, char* weights);
processor_graph_access create3DGrid (int size, char* weights);
void allPathsShortestPaths (vector< vector<int> >* predecessor_matrix, vector< vector<int> >* distance_matrix, graph_access* G);

#endif /* PROCESSOR_GRAPH_H */
