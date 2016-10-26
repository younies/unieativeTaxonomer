//
//  testingYRJ2.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/21/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "testingYRJ.hpp"


void Tester::testingGenomeLevelWithNewMethodology(YRJObject * yrj   , int differences)
{
    //using my new methodology to make it better
    
    auto  unieativeHits = this->getUnieativeHitsGenus(yrj ,  differences );
    
    auto finalUnieative = this->unieativeLCAKraken(unieativeHits);
    
    auto level = this->getLeastCommonLevel(yrj, finalUnieative);
    
    
    if(this->finalResult.count(level))
        finalResult[level] ++;
    else
        finalResult[level] = 1;
    
}



void Tester::testingGenomeLevel(YRJObject * yrj   , int differences)
{
    auto hitNumbers = this->getKrakenLCAs(yrj, differences);
    
    //to get kraken Final Taxonomy
    
    if(hitNumbers.size() == 0)
    {
        this->finalResult[this->notConsidered] ++;
        return;
    }
    
    auto krakenShort = this->pruinedTree->getTheMaximumKRAKENhit(hitNumbers);
    
    auto level = this->getLeastCommonLevel(yrj, krakenShort);
    
    if(this->finalResult.count(level))
        finalResult[level] ++;
    else
        finalResult[level] = 1;
    
}




bool Tester::isKrakenCatch(YRJObject * yrj )
{
    //to find all the LCAs
    for(auto kmer: yrj->kmersVector)
    {
        auto hits = this->hits_kmer_with_differences(kmer, 0);
        
        if(hits.size() != 0)
            return true;
    }
    
    return false;

}

map<short, int>  Tester::getKrakenLCAs(YRJObject * yrj , int differences)
{
    map<short, int> hitNumbers;
    //to find all the LCAs
    for(auto kmer: yrj->kmersVector)
    {
        auto hits = this->hits_kmer_with_differences(kmer, differences);
        
        if(hits.size() == 0){
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
    return hitNumbers;
}



map<short, int>  Tester::getUnieativeHitsGenus(YRJObject * yrj , int differences )
{
    map<short, int> unieativeHits;
    
    
    
    for(auto kmer : yrj->kmersVector)
    {
        auto hits = this->hits_kmer_with_differences(kmer, differences);
        
        if(hits.size() == 0){
            continue;
        }
        
        set<short> genusHits;
        
        for (auto hit: hits)
        {
            genusHits.insert( this->pruinedTree->getGenusParent(hit));
        }
        
        
        for(auto genusHit : genusHits)
        {
            if(unieativeHits.count(genusHit))
                unieativeHits[genusHit]++;
            else
                unieativeHits[genusHit] = 1;
        }
        
    }
    
    return unieativeHits;
     
}



short Tester::unieativeLCAKraken(map<short, int> unieativeHits)
{
    int maximum = 0;
    
    vector<short> maxNodes;
    for(auto hit: unieativeHits)
        maximum = max(maximum , hit.second);
    
    for(auto unieative: unieativeHits)
        if(unieative.second == maximum)
            maxNodes.emplace_back(unieative.first);
    
    
    return this->pruinedTree->getGlobalLCA(maxNodes);

    
}



string Tester::getLeastCommonLevel(YRJObject * yrj, short krakenShort)
{
    auto krakenUID   = this->pruinedTree->getTheUIDFromShort(krakenShort);
    
    auto krakenNode  = this->bigTree->getNodeFromIndex(this->bigTree->uid_to_index(krakenUID));
    
    auto indexYRJ = this->bigTree->uid_to_index(yrj->uid);
    
    auto yrjNode     = this->bigTree->getNodeFromIndex(indexYRJ);
    
    auto testingLevelIndex = this->bigTree->get_LCA_between_Two_Nodes(krakenNode, yrjNode);
    
    auto levelUID = this->bigTree->getNodeFromIndex(testingLevelIndex);
    
    string level = this->bigTree->get_level(levelUID);
    
    return level;

}


