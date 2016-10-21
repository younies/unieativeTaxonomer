//
//  testingYRJ2.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/21/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "testingYRJ.hpp"



void Tester::testingGenomeLevel(YRJObject * yrj   , int differences)
{
    map<int, int> hitNumbers;
    
    for(auto kmer: yrj->kmersVector)
    {
        auto hits = this->hits_kmer_with_differences(kmer, differences);
        
        auto lca = this->pruinedTree->getGlobalLCA(hits);
        
        if(hitNumbers.count(lca))
            hitNumbers[lca]++;
        else
            hitNumbers[lca] = 1;
    }
    
    
}
