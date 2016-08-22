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
#include "helpers.hpp"
#include "Tree.hpp"

class CoreTaxonomer
{
    LONGS coreHashNodesSize = 0;
    vector<HashedNode> coreHashedNodes;
    vector<YRJObject *> yrjVector;
    
    LONGS kmerLength = 31; // need some changes
    
    LONGS startIndex = 0;
    
    bitset<64> hash;
    int sizeOfHash = 64;
    int sizeOfHahsedRawPart = 32;
    LONG reverseHashBits;
    //void copyYRJUnieativeInside(YRJUnieative &yrjUnieative);
    
    Tree * globalTree;
    
public:
    CoreTaxonomer( vector<YRJObject *> yrjVector , string hash);
    ~CoreTaxonomer();
    
    pair<LONGS, LONGS> getThePlaceOfKmer(pair<SHORT, SHORT>  rawKmer);
    void fillAllTheCoreData();
    
    HashedNode getTheHashedKmer(LONG kmer);// return the rawKmer and the hashedKmer from a kmer
    
    void updateHashValue(string hash);
    
    LONG reverseKmer(LONG kmer);
    
    vector<pair< short , short> > getShortNameFromKmer(LONG kmer); // return the short names associated with the number of differences
    
    pair<short, short>  scanAtIndex( LONGS index , pair<SHORT, SHORT> pairHashed);
    
    pair<SHORT , SHORT> convertINTtoPairShort(INT kmerINT);
    
    void writeTheCoreData(string path);
    
};

#endif /* CorseTaxonomer_hpp */
