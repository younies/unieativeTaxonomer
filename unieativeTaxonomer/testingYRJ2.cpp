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
    map<short, int> hitNumbers;
    
    
    //to find all the LCAs
    for(auto kmer: yrj->kmersVector)
    {
        auto hits = this->hits_kmer_with_differences(kmer, differences);
        
        if(hits.size() == 0){
            cout << "hits 0\n";
            continue;
        }
        if(hits.size() > 0)
        {
            auto lca = this->pruinedTree->getGlobalLCA(hits);
        
            if(hitNumbers.count(lca))
                hitNumbers[lca]++;
            else
                hitNumbers[lca] = 1;
        }
        
        cout << "passed hash\n";
    }
    
    
    cout << "passed 1 \n";
    
    
    //to get kraken
    
    cout << hitNumbers.size() << endl;
    auto krakenShort = this->pruinedTree->getTheMaximumKRAKENhit(hitNumbers);
    
    cout << "passed 2 \n";

    
    auto krakenUID   = this->pruinedTree->getTheUIDFromShort(krakenShort);
    
    cout << "passed 3 \n";

    auto krakenNode  = this->bigTree->getNodeFromIndex(this->bigTree->uid_to_index(krakenUID));
    cout << "passed 4 \n";

    cout << yrj->uid << endl;
    
    auto index = this->bigTree->uid_to_index(yrj->uid);
   cout << "passed 4.1 \n";
    auto yrjNode     = this->bigTree->getNodeFromIndex(index);
    
    cout << "passed 5 \n";

    
    auto testingLeveLUID = this->bigTree->get_LCA_between_Two_Nodes(krakenNode, yrjNode);
    
    
    auto levelUID = this->bigTree->getNodeFromIndex(this->bigTree->uid_to_index(testingLeveLUID));
    cout << testingLeveLUID << "   " << this->bigTree->get_level(levelUID);
    
    
    
    
}
