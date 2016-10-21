//
//  tree_support.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/21/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "Tree.hpp"



//implementation for maximum hits
short Tree::getTheMaximumKRAKENhit(map<short, int> & originalHitsMap)
{
    auto retHitsMap = originalHitsMap;
    
    
    //for aggregating hits
    for(auto hit : originalHitsMap)
    {
        auto lca = hit.first;
        auto parent = hit.first;
        do
        {
            parent = this->getParentShortName(parent);
            if(originalHitsMap.count(parent))
                retHitsMap[lca] += originalHitsMap[parent];
            
        }while(parent > 0);
    }

    
    //for finding the maximum hits
    int maxHits = 0;
    for(auto hit : retHitsMap)
        if(hit.second > maxHits)
            maxHits = hit.second;
    
    
    //for extracting the maximum nodes
    vector<short> ret;
    for(auto hit : retHitsMap)
        if(hit.second == maxHits)
            ret.emplace_back(hit.first);
    
    return this->getGlobalLCA(ret);
    
}
