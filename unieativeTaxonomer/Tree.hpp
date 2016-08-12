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

class Tree {
    vector< TreeNode * > treeNodesVector;
    vector<LONGS> fromShortNameToUid;
    
public:
    Tree(string pathToTheTree);
    
    void setChildrenForTheNodes();
    
};
#endif /* Tree_hpp */
