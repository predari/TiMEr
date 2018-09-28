/* Roland Glantz <roland.glantz@kit.edu> */

#ifndef ALMA_GRAPH_H
#define ALMA_GRAPH_H

using namespace std;
#include "../lib/data_structure/graph_access.h"
#include "./alma/permutation/permutation.h"
#include <stdlib.h>
#include <string>
#include <vector>
#include <unordered_map>

class alma_graph_access : public graph_access
{
protected:
int64_t numL; ///< number of labels >= number of vertices
int64_t numDigitsL; ///< number of digits in binary representation of labels
int64_t numDigitsBlockID; ///< number of digits for block ID
int64_t numDigitsForID; ///< number of digits for identifying vertices within a block
vector<int64_t> vlabels; ///<vertex labels
unordered_map<int64_t, NodeID> label2vertex;

public:

alma_graph_access();

~alma_graph_access();

/**
 * @return numDigitsL.
 */
inline int64_t getNumDigitsL() {
return numDigitsL;
}

/**
 * @return numDigitsBlockID.
 */
inline int64_t getNumDigitsBlockID() const {
return numDigitsBlockID;
}

/**
 * @return numDigitsForID (numDigitsL = numDigitsBlockID + numDigitsForID).
 */
inline int64_t getNumDigitsForID() const {
return numDigitsForID;
}

/**
 * setting numDigitsL.
 */
void setNumDigitsL(int64_t ndl);
    
/**
 * @return numL.
 */
inline int64_t getNumL() const {
return numL;
}
    
/**
 * setting numL.
 */
void setNumL(int64_t nl);

/**
 * Returns label of vertex
 * @param v vertex.
 * @return label of vertex v
 */
inline int64_t getLabel(int64_t v) {
return(vlabels[v]);
}

/**
 * Assigns label l to vertex v
 * @param[in] vertex v.
 * @param[out]label l
 */
inline void setLabel(NodeID v, int64_t l) {
vlabels[v] = l;
}

/**
 * Resizes vlabels
 * @param[in] size.
 */
inline void resize_vlabels(NodeID size) {
vlabels.resize(size);
}

/**
 * Finds initial labels for the vertices based on a given partition
 * param[out] numDigitsBlockID
 * param[out] numDigitsForID
 * param[out] numDigitsL
 * param[out] numL
 */
void getVertexLabelsFromPartition(alma_graph_access & g);

/**
 * Permutes vertex labels according to permutation pi.
 * @param[in] permutation pi.
 */
void permuteVertexLabels(alma_graph_access & g, Permutation pi);

/**
 * Builds vector label2vertex
 */
void makeLabel2Vertex(alma_graph_access & g);
    
/**
 * Writes vertices, vertex labels and edges to standard out
 */
void display(alma_graph_access & g);
    
/**
 * For all vertices in graph, calculate distance to vertex v.
 * @param[in] vertex v.
 * @param[in] vector of distances to v, dist2v.
 * @param[out] dist2v.
 */
void getDistances2vertex(alma_graph_access & g, vector<int64_t> &dist2v, int64_t v);
    
/**
 * @returns true if the graph is a partial cube and false, otherwise
 * if graph is a partial cube, assigns Hamming labels to vertices (in vlabels[])
 */
bool getLabelsPartialCube(alma_graph_access & g);
    
/**
 * Writes number of missfits at positions of labels to standard out
 */
void statisticsOnMissfits(alma_graph_access & g);

/**
 * Returns vertex with label l
 * @param[in] l label
 * @return vertex with label l (-1 if there is no such vertex)
 */
NodeID getVertexWithLabel(int64_t l);
    
/**
 * Get total cut of graph, i.e., sum of edge weights beween parts.
 * @return total cut.
 */
int64_t getTotalCut(alma_graph_access & g);
    
/**
* Check whether partition, as given by vlabels, is valid.
* @return true iff partition is valid.
*/
bool checkPartition(alma_graph_access & g, int64_t numB);
    
/**
 * Get total communcation cost, i.e., sum of Hamming distances between edges' end vertices.
 * @return total communication cost.
 */
int64_t getTotalCommunicationCost(alma_graph_access & g);


double  getMeanDegree(alma_graph_access & g);
int64_t getMaxDegree(alma_graph_access & g);
int64_t getMinDegree(alma_graph_access & g);
  
int64_t getMaxCommunicationCost(alma_graph_access & g);
  
/**
 * Get total communcation cost, i.e., sum of Hamming distances between edges' end vertices.
 * @return total communication cost.
 */

  
int64_t getTotalHammingDistance(alma_graph_access & g);


/**
 * Contracts graph g.
 * @param[in] g
 * @param[out] parent
 * parent[v] = parent of v in contracted graph for all vertices v of g
 * @return contracted g
 */
 alma_graph_access* contract(alma_graph_access* h_ptr, alma_graph_access & g, vector<NodeID> & parents);

/**
 * Checks neighborhood of v for gains due to
 * label change from even to odd.
 * @param[in] vertex v
 * @return potential gain around v
 */
int64_t gainSwapUp(alma_graph_access & g, int64_t v);
    
/**
 * Checks neighborhood of w for gains due to
 * label change from odd to even.
 * @param[in] vertex w
 * @return potential gain around w
 */
int64_t gainSwapDown(alma_graph_access & g, int64_t w);

/**
 * Get weight of edge between v and w, if any. Zero otherwise.
 * @param[in] vertex v
 * @param[in] vertex w
 * @return weight of edge {v, w} if this edge exists, 0 otherwise
 */
EdgeWeight getWeightBetween(alma_graph_access & g, int64_t v, int64_t w);
    
/**
 * Swaps vertex labels for intensification gains.
 */
 void swapVertexLabelsIntensification(alma_graph_access & g, int64_t digit);

/**
 * Swaps vertex labels for diversification gains.
 */
void swapVertexLabelsDiversification(alma_graph_access & g, int64_t digit);

};

#endif /* ALMA_GRAPH_H */
