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
    
    cout << "passed 111 \n";

    //to find all the LCAs
    for(auto kmer: yrj->kmersVector)
    {
        cout << "passed 112 \n";

        auto hits = this->hits_kmer_with_differences(kmer, differences);
        
        cout << "passed 113 \n";

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
    
    cout << "passed 22 \n";

    
    auto krakenUID   = this->pruinedTree->getTheUIDFromShort(krakenShort);
    cout << "passed 23 \n";
    cout << krakenUID << endl;


    auto krakenNode  = this->bigTree->getNodeFromIndex(this->bigTree->uid_to_index(krakenUID));

    cout << "passed 24 \n";

    cout << yrj->uid << endl;
    
    auto indexYRJ = this->bigTree->uid_to_index(yrj->uid);
   
    
     cout << indexYRJ << endl;
    cout << "passed 25 \n";

    auto yrjNode     = this->bigTree->getNodeFromIndex(indexYRJ);
    
    cout << "passed 26 \n";

    
    auto testingLevelIndex = this->bigTree->get_LCA_between_Two_Nodes(krakenNode, yrjNode);
    cout << "passed 27 \n";
    
    
    cout <<  testingLevelIndex << endl;
    
    
    
    auto levelUID = this->bigTree->getNodeFromIndex(testingLevelIndex);
    cout << this->bigTree->get_level(levelUID);
    
    
    
    
}
