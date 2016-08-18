//
//  HashedNode.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef HashedNode_hpp
#define HashedNode_hpp

#include "headers.h"

//abstracting the new hashed kmer
struct HashedNode{
    short index;
    pair<SHORT, SHORT> rawKmer; //the non hashed part in the kmer
    pair<SHORT, SHORT> hashedKmer ; // the part that is hashed in the kmer
    
    HashedNode(){
        index = 0;
        rawKmer.first = 0;
        rawKmer.second = 0;
        hashedKmer.first = 0;
        hashedKmer.second = 0;
    }
};

bool  hashedNodeCompare( const HashedNode &lhs, const HashedNode &rhs);
//end of the hasehd kmer implementation

#endif /* HashedNode_hpp */
