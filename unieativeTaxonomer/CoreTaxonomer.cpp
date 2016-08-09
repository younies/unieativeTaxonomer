//
//  CorseTaxonomer.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "CoreTaxonomer.hpp"


CoreTaxonomer::CoreTaxonomer(vector<YRJUnieative> yrjUnieativeVector)
{
    //to calculate the size of the whole database hashed kmers
    this->coreHashNodesSize = 0;
    for(YRJUnieative node : yrjUnieativeVector)
        this->coreHashNodesSize += node.hashedKmersSize;
    
    //to build the hashed database
    this->coreHashedNodes.resize(this->coreHashNodesSize);
}

CoreTaxonomer::~CoreTaxonomer()
{
    
    
}