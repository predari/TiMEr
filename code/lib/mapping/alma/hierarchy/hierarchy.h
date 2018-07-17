/*
 * hierarchy.h
 *
 *  Created on: Dec 22, 2016
 *      Author: Roland
 */

#ifndef HIERARCHY_H
#define	HIERARCHY_H

#include "../../alma_graph_access.h"
#include <iostream>
#include <vector>

using namespace std;

class Hierarchy {
private:
    int64_t numLevels;
    vector<alma_graph_access*> graphsVec;
    vector<vector<NodeID>> parentsVec;
    
public:

    Hierarchy(int64_t numLevels);

    ~Hierarchy(void);
    
    /**
     * Get number of levels in the hierarchy
     * @return numLevels.
     */
    int64_t getNumberOfLevels(void);

    /**
     * Get number of vertices in graph on level l.
     * param[in] level l.
     * @return number of vertices in graph on level l.
     */
    NodeID getNumberVerticesOnLevel(int64_t l);
    
    /**
     * display parents on level l (standard out).
     * param[in] level l.
     */
    void displayParentsOnLevel(int64_t l);
    
    /**
     * Insertion of newG at position pos of hierarchy.
     * @param[in] position pos.
     */
    void insertGraph(alma_graph_access* newG, int64_t pos);

    /**
     * Return finest graph of hierarchy.
     * @return graphsVec[0].
     */
    alma_graph_access retFinestGraph(void);

    /**
     * Insertion of parents at position pos of hierarchy.
     * @param[in] position pos.
     */
    void insertParents(vector<NodeID>& parents, int64_t pos);
    
    /**
     * Update vertex labels of finest graph in hierarchy.
     * Graph on finest level is modfied
     */
    void updateLabelsOnFinestLevel(void);
};

#endif	/* HIERARCHY_H */

