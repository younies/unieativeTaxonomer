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

        cout << kmer << endl;
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
        
    }
    
    
    
    
    //to get kraken
    
    if(hitNumbers.size() == 0)
    {
        cout << "not_considered\n";
        return;
    }
    auto krakenShort = this->pruinedTree->getTheMaximumKRAKENhit(hitNumbers);
    

    
    auto krakenUID   = this->pruinedTree->getTheUIDFromShort(krakenShort);
    cout << krakenUID << endl;


    auto krakenNode  = this->bigTree->getNodeFromIndex(this->bigTree->uid_to_index(krakenUID));


    cout << yrj->uid << endl;
    
    auto indexYRJ = this->bigTree->uid_to_index(yrj->uid);
   
    
     cout << indexYRJ << endl;

    auto yrjNode     = this->bigTree->getNodeFromIndex(indexYRJ);
    

    
    auto testingLevelIndex = this->bigTree->get_LCA_between_Two_Nodes(krakenNode, yrjNode);
    
    
    cout <<  testingLevelIndex << endl;
    
    
    
    auto levelUID = this->bigTree->getNodeFromIndex(testingLevelIndex);
    cout << this->bigTree->get_level(levelUID) << endl;
    
    
    
    
}
