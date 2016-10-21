//
//  Tree.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/11/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef Tree_hpp
#define Tree_hpp
#include "YRJObject.hpp"
#include "headers.h"
#include "TreeNode.hpp"

class Tree
{
    vector< TreeNode * > treeNodesVector;
    //vector<LONGS> fromShortNameToUid;
    LONGS treeSize = 11000;
    vector<LONGS> numberOfLeaves;
    void connectChildren();
    LONGS setNumberOfLeaves(TreeNode * node);
    
    

public:
    short root = 0; /// should always hold the short for the root
    Tree(string pathToTheTree);
    
    short getParentShortName(short shortName);
    
    void buildTheNumberOfLeaves( );
    
    //to get the g_statistics domenrator
    
    
    short getTowLCA(short first , short second);
    
    short getGlobalLCA( vector<short> &  hitted_nodes);
    
    short getGlobalLCA( vector<pair< short , short> > & vectorOfResults);// take vector of pairs of indices and differences

    LONGS getNumberOfLeaves(short shortName);
    
    short getLCA_for_child_node( TreeNode * node , vector<short> & hitted_nodes );
    
    G_Statistics calculateG_Statistics(LONG kmer , short nodeShortName  , vector<short> & hitted_nodes );
    
    LONGS getTheUIDFromShort(short shortName);
    
    vector<YRJObject *> getYRJobjects(string path );
    
    
    
    short getTheMaximumKRAKENhit(map<short, int> & hitNumbers);
    
};
#endif /* Tree_hpp */
