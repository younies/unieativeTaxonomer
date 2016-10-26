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
    
    short genus; // the short id for genus level of this node
    short species; // the short id for the genus level of this node
    TreeNode * parent;
    vector<TreeNode* > children;
    short shortName;
    short parentShortName;
    bool tagged;
    LONGS uid, parentUid;
    
    
    short speciesParent = -2;
    short genusParent = -2;
    

    
    

    TreeNode(LONGS uid , LONGS parentUid , short shortName , short parentShortName , bool tagged);
};

#endif /* TreeNode_hpp */
