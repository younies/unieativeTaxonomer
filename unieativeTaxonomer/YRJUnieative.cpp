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
    vector<char> newHash;
    
    for(int i = 0 , n = (int)hash.size() ; i < n ; ++i)
    {
        newHash.push_back(hash[i]);
        newHash.push_back(hash[i]);
    }
    string newStringHash( newHash.begin() , newHash.end());
    this->hash = newStringHash;
}


//building the converter




HashedNode YRJUnieative::convetToHashed(LONG kmer)
{
    HashedNode node;
    
    node.index = this->getIndex();
    
    
    return node;
}