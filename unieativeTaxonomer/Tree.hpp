//
//  Tree.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/11/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef Tree_hpp
#define Tree_hpp

#include "headers.h"
#include "TreeNode.hpp"
#include "YRJUnieative.hpp"

class Tree
{
    vector< TreeNode * > treeNodesVector;
    vector<LONGS> fromShortNameToUid;
    LONGS treeSize = 11000;
    vector<LONGS> numberOfLeaves;
    void connectChildren();
    LONGS setNumberOfLeaves(TreeNode * node);

public:
    Tree(string pathToTheTree);
    
    short getParentShortName(short shortName);
    
    void buildTheNumberOfLeaves( );
    
    //to get the g_statistics domenrator
    
    short getTowLCA(short first , short second);
    
    short getGlobalLCA( vector<pair< short , short> > & vectorOfResults);// take vector of pairs of indices and differences

    LONGS getNumberOfLeaves(short shortName);
    
    short getLCA_for_child_node( TreeNode * node , vector<short> & hitted_nodes );
    
    G_Statistics calculateG_Statistics(LONG kmer , short nodeShortName  , vector<short> & hitted_nodes );
    
    vector<YRJUnieative *> getYRJUnieariveVector(string path ,string hash);
};
#endif /* Tree_hpp */
