//
//  TreeNode.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/11/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//


/**
 This class represinting a single node
 from the pruning phylogenetic tree.
 */

#ifndef TreeNode_hpp
#define TreeNode_hpp

#include "headers.h"


class TreeNode
{

public:
    TreeNode * parent;
    list<TreeNode* > children;
    short shortName;
    bool tagged;
    LONG uid, parentUid;
    
    

    TreeNode(LONG uid , LONG parentUid , short shortName , bool tagged);
};

#endif /* TreeNode_hpp */
