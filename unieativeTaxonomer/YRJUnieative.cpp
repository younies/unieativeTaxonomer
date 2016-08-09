//
//  YRJUnieative.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/9/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "YRJUnieative.hpp"



YRJUnieative::YRJUnieative(string path , string hash): YRJObject(path)
{
    this->hash = hash;
}


//building the converter


unsigned int 

HashedNode YRJUnieative::convetToHashed(LONG kmer)
{
    HashedNode node;
    
    node.index = this->getIndex();
    
    
    
    
    return node;
}