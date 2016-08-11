//
//  CorseTaxonomer.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef CoreTaxonomer_hpp
#define CoreTaxonomer_hpp
#include "HashedNode.hpp"
#include "headers.h"
#include "YRJUnieative.hpp"

class CoreTaxonomer
{
    LONGS coreHashNodesSize = 0;
    vector<HashedNode> coreHashedNodes;
    vector<YRJUnieative *> yrjUnieativeVector;
    LONGS startIndex = 0;
    void copyYRJUnieativeInside(YRJUnieative &yrjUnieative);
    
public:
    CoreTaxonomer( vector<YRJUnieative *>  &yrjUnieativeVector);
    ~CoreTaxonomer();
    
    pair<LONGS, LONGS> getThePlaceOfKmer(INT rawKmer);
    void fillAllTheCoreData();
    
    pair<INT, INT> getTheHashedKmer(LONG kmer);// return the rawKmer and the hashedKmer from a kmer
    
    
    
};

#endif /* CorseTaxonomer_hpp */
