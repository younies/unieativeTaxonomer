//
//  HashedNode.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "HashedNode.hpp"

bool hashedNodeCompare(  const HashedNode &lhs, const HashedNode &rhs)
{
    return std::tie(lhs.rawKmer, lhs.index , lhs.hashedKmer ) < std::tie(rhs.rawKmer, rhs.index , rhs.hashedKmer);
} // for the binartysearch or the sorting algorithms


