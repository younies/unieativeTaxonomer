//
//  TestingKrakenInDifferenceLevels.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/9/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "TestingKrakenInDifferenceLevels.hpp"



void testOutYRJfile( ofstream * writeFile , vector<Hash*> &hashes, Tree * tree, BigTree * bigTree, YRJObject * yrjObject  ,  Unieative * unieative ,int numOfDifferences )
{
    
    
    for (int  currDiff = 0 ; currDiff <= numOfDifferences ; ++currDiff) {
        LONGS assignedUid = unieative->getFinalUIDs(yrjObject, currDiff);
        
        if(bigTree->isBothInGenusLevel(assignedUid, <#LONGS uid_second#>))
    }
    
    
    
}
