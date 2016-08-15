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
#include "helpers.hpp"
#include "Tree.hpp"

class CoreTaxonomer
{
    LONGS coreHashNodesSize = 0;
    vector<HashedNode> coreHashedNodes;
    vector<YRJUnieative *> yrjUnieativeVector;
    
    LONGS startIndex = 0;
    string hash;
    LONG reverseHashBits;
    void copyYRJUnieativeInside(YRJUnieative &yrjUnieative);
    
    Tree * globalTree;
    
public:
    CoreTaxonomer( vector<YRJUnieative *>  &yrjUnieativeVector , string hash);
    ~CoreTaxonomer();
    
    pair<LONGS, LONGS> getThePlaceOfKmer(INT rawKmer);
    void fillAllTheCoreData();
    
    pair<INT, INT> getTheHashedKmer(LONG kmer);// return the rawKmer and the hashedKmer from a kmer
    
    void updateHashValue(string hash);
    
    LONG reverseKmer(LONG kmer);
    
    vector<pair< short , short> > getShortNameFromKmer(LONG kmer); // return the short names associated with the number of differences
    
    pair<short, short> scanAtIndex(LONGS index ,INT hashed); //returns the short name and
    
    short getTowLCA(short first , short second);
    
    short getGlobalLCA( vector<pair< short , short> > vectorOfResults);// take vector of pairs of indices and differences
    
};

#endif /* CorseTaxonomer_hpp */
