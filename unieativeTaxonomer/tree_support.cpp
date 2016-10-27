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







short Tree::getSpeciesParent(YRJObject * yrj)
{
    if(yrj->speciesParent != -2)
        return yrj->speciesParent;
    
    
    
    
    auto speciesParent = yrj->shorName;
    
    while (speciesParent != 0)
    {
        
        auto uid = this->treeNodesVector[speciesParent]->uid;
        auto index_in_BigTree = this->bigTree->uid_to_index(uid);
        auto level = this->bigTree->getLevel(index_in_BigTree);
        
        if(level == this->species)
        {
            yrj->speciesParent = speciesParent;
            return speciesParent;
        }
        
        speciesParent = this->treeNodesVector[speciesParent]->parentShortName;
    }
    
    yrj->speciesParent = -1;
    return -1;
    
    
}



short  Tree::getGenusParent(YRJObject * yrj)
{
    cout << "genusParent" << endl;
    if(yrj->genusParent != -2)
        return yrj->genusParent;
    
    
    auto genusParent = yrj->shorName;
    
    cout << genusParent << endl;
    
    while (genusParent != 0)
    {
        
        cout << genusParent << endl;
        
        auto uid = this->treeNodesVector[genusParent]->uid;
        auto index_in_BigTree = this->bigTree->uid_to_index(uid);
        auto level = this->bigTree->getLevel(index_in_BigTree);
        
        if(level == this->genus)
        {
            yrj->genusParent = genusParent;
            return genusParent;
        }
        
        genusParent = this->treeNodesVector[genusParent]->parentShortName;

    }
    
    yrj->genusParent = -1;
    return -1;
}


short Tree::getSpeciesParent(short shortName)
{
    if( this->treeNodesVector[shortName]->speciesParent != -2)
        return this->treeNodesVector[shortName]->speciesParent;
    
    
    
    
    auto speciesParent = shortName;
    
    while (speciesParent != 0)
    {
        auto uid = this->treeNodesVector[speciesParent]->uid;
        auto index_in_BigTree = this->bigTree->uid_to_index(uid);
        auto level = this->bigTree->getLevel(index_in_BigTree);
        
        if(level == this->species)
        {
            this->treeNodesVector[shortName]->speciesParent = speciesParent;
            return speciesParent;
        }
        
        speciesParent = this->treeNodesVector[speciesParent]->parentShortName;
    }
    
    this->treeNodesVector[shortName]->speciesParent = -1;
    return -1;
    
}


short Tree::getGenusParent(short shortName)
{
    if(this->treeNodesVector[shortName]->genusParent != -2)
        return this->treeNodesVector[shortName]->genusParent;
    
    
    auto genusParent = shortName;
    
    while (genusParent != 0)
    {
        auto uid = this->treeNodesVector[genusParent]->uid;
        auto index_in_BigTree = this->bigTree->uid_to_index(uid);
        auto level = this->bigTree->getLevel(index_in_BigTree);
        
        if(level == this->genus)
        {
            this->treeNodesVector[shortName]->genusParent = genusParent;
            return genusParent;
        }
        
        genusParent = this->treeNodesVector[genusParent]->parentShortName;
        
    }
    
    this->treeNodesVector[shortName]->genusParent = -1;
    return -1;
}

