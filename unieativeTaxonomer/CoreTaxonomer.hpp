//
//  CorseTaxonomer.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef CoreTaxonomer_hpp
#define CoreTaxonomer_hpp
#include "HashedNode.hpp"
#include "headers.h"
#include "YRJUnieative.hpp"

class CoreTaxonomer
{
    LONGS coreHashNodesSize = 0;
    vector<HashedNode> coreHashedNodes;
    
    
    
public:
    CoreTaxonomer( vector<YRJUnieative> yrjUnieativeVector);
    ~CoreTaxonomer();
    
    
};

#endif /* CorseTaxonomer_hpp */
