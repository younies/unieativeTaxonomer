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
    
};
#endif /* Tree_hpp */
