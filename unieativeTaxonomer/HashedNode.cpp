//
//  HashedNode.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "HashedNode.hpp"

bool  * hashedNodeCompare( HashedNode &lhs, const HashedNode &rhs) {
    
    bool ret = (lhs < rhs);
    return &ret; } // for the binartysearch or the sorting algorithms
