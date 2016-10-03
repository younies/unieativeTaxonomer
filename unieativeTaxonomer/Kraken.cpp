//
//  Kraken.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/3/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "Kraken.hpp"



SHORT Unieative::getLCA(LONG kmer , int & differences)
{
    
    
    unordered_set<short> hits;
    for(auto hash: this->Hashes)
    {
        auto tempHits = getNumberOfDifference(hash->getIndexStream(), hash->getDataStream(), hash, kmer);
        
        for(auto pair : tempHits)
            if(pair.second <= differences)
                hits.insert(pair.first);
        
    }
    
    vector<short> uniqueHits(hits.size());
    
    copy( hits.begin() , hits.end() , uniqueHits.begin() );
    
    return this->tree->getGlobalLCA(uniqueHits);
    
    
}








vector<INT> Unieative::getFinalUIDs(YRJObject * yrjObject , int differnces)
{
    unordered_map<SHORT , int> LCAs;
    
    //find the number of hits for each lca corresponding to each kmer
    for(auto kmer : yrjObject->kmersVector)
    {
        auto  lca  = this->getLCA(kmer, differnces);
        LCAs.count(lca) ? ++LCAs[lca] : LCAs[lca] = 1;
    }
    
    
    auto comulatedLCAS = LCAs;
    // find the total number of hits in the parents too
    for(auto child : LCAs)
    {
        short parent = tree->getParentShortName( child.first);
        while ( parent != tree->root){
            if(LCAs.count(parent))
               comulatedLCAS[child.first] += LCAs[parent];// add to the child more hits
            parent = tree->getParentShortName(parent);
        }
    }
    
    vector<INT> ret;
    if(LCAs.empty())
        return ret;
    
    int maxi = -1;
    for(auto cLCA: comulatedLCAS )
        if(cLCA.second > maxi)
            maxi = cLCA.second;
    
    for(auto cLCA : comulatedLCAS)
        if(cLCA.second == maxi)
            ret.emplace_back(cLCA.first);
    
    return ret;
}











