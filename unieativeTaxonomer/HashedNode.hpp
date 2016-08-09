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
    unsigned int rawKmer; //the non hashed part in the kmer
    unsigned int hashedKmer; // the part that is hashed in the kmer
    
    bool operator< ( const HashedNode& y) {
        return std::tie(this->rawKmer, this->index , this->hashedKmer ) < std::tie(y.rawKmer, y.index , y.hashedKmer);
    }
};

bool hashedNodeCompare(HashedNode &lhs, HashedNode &rhs);
//end of the hasehd kmer implementation

#endif /* HashedNode_hpp */