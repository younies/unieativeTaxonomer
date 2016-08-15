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
    
public:
    Tree(string pathToTheTree);
    
    
    short getParentShortName(short shortName);
    
};
#endif /* Tree_hpp */
