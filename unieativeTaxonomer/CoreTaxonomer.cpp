//
//  CorseTaxonomer.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "CoreTaxonomer.hpp"


CoreTaxonomer::CoreTaxonomer(vector<YRJUnieative *> & yrjUnieativeVector)
{
    this->yrjUnieativeVector = yrjUnieativeVector;
    //to calculate the size of the whole database hashed kmers
    this->coreHashNodesSize = 0;
    for(YRJUnieative* node : yrjUnieativeVector)
        this->coreHashNodesSize += node->hashedKmersSize;
    
    //to build the hashed database
    this->coreHashedNodes.resize(this->coreHashNodesSize);
    
    //compy files in the coreHahsdNodes
    this->fillAllTheCoreData();
    
    
    //sort all the core data
    sort(this->coreHashedNodes.begin(), this->coreHashedNodes.end());
    
}



CoreTaxonomer::~CoreTaxonomer()
{
    
    
}




//merge yrjUnieative inside the the core data

void CoreTaxonomer::copyYRJUnieativeInside(YRJUnieative & yrjUnieative)
{
    yrjUnieative.fillTheHashedNodesVector();
    for(LONGS i = 0  , n = yrjUnieative.hashedKmers.size() ; i < n ; ++ i)
        this->coreHashedNodes[this->startIndex++] = yrjUnieative.hashedKmers[i];
    yrjUnieative.clearAllTheData();
}



void CoreTaxonomer::fillAllTheCoreData()
{
    this->startIndex = 0;
    for (LONGS i = 0 , n = this->yrjUnieativeVector.size(); i < n ; ++i)
    {
        this->copyYRJUnieativeInside(*(this->yrjUnieativeVector[i]));
    }
    
    this->startIndex = 0;
    
}


pair<LONGS, LONGS>  CoreTaxonomer::getThePlaceOfKmer(INT rawKmer)
{
    LONGS start = 0 , end = this->coreHashNodesSize - 1 , mid  = (start + end)/2;
    
    while (end >= start)
    {
        mid  = (start + end)/2;
        if(this->coreHashedNodes[mid].rawKmer == rawKmer)
            break;
        else if (this->coreHashedNodes[mid].rawKmer > rawKmer)
            end = mid - 1;
        else
            start = mid + 1;
    }
    
    
    LONGS newEndStart =  mid ,  newEndEnd  = mid;
    
    //for setting the start
    while ( newEndStart > start )
    {
        mid = (start + newEndStart) / 2;
        
        if(this->coreHashedNodes[mid].rawKmer == rawKmer)
            newEndStart = mid;
        else
            start = mid + 1;
    }
    
    //for setting the end
    
    while (end > newEndEnd)
    {
        mid = (end + newEndEnd)/2;
        if(this->coreHashedNodes[mid].rawKmer == rawKmer)
            newEndEnd = mid;
        else
            end = mid  - 1;
        
    }
    
    return make_pair(start, end);
    
}



