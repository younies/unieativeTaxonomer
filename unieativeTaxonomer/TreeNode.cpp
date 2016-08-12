//
//  TreeNode.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/11/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "TreeNode.hpp"


TreeNode::TreeNode(LONG uid , LONG parentUid , short shortName , bool tagged)
{
    this->uid           = uid;
    this-> parentUid    = parentUid;
    this->shortName     = shortName;
    this->tagged        = tagged;
}

