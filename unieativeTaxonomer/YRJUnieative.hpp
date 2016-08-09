//
//  YRJUnieative.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/9/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//


/**
 
 This class is for having the hashed kmer as that specified in the documents
 
 */

#ifndef YRJUnieative_hpp
#define YRJUnieative_hpp

#include "headers.h"
#include "YRJObject.hpp"
#include "HashedNode.hpp"


class YRJUnieative: YRJObject
{
    vector< HashedNode > hashedKmers;
    string hash;
    HashedNode convetToHashed(LONG kmer);
    void extractKeys( HashedNode &node , LONG & kmer);
public:
    YRJUnieative(string path , string hash);
    ~YRJUnieative();
    
    LONGS hashedKmersSize;
    
    void fillTheHashedNodesVector();

};


#endif /* YRJUnieative_hpp */
