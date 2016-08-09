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
    //to find the new hash
    vector<char> newHash;
    for(int i = 0 , n = (int)hash.size() ; i < n ; ++i)
    {
        newHash.push_back(hash[i]);
        newHash.push_back(hash[i]);
    }
    string newStringHash( newHash.begin() , newHash.end());
    this->hash = newStringHash;
    
}



YRJUnieative::~YRJUnieative()
{
    this->hashedKmers.clear();
    this->hashedKmersSize = 0;
    cout << "yrj unieative is destroyed\n";
}

//building the converter

void YRJUnieative::extractKeys( HashedNode &node , LONG & kmer)
{
    unsigned int idKmer = 0 , hashed = 0;
    int n = (int)this->hash.size();
    unsigned long one = 1;
    for (int i = 0 ; i < n ; ++i)
    {
        if(hash[i] == '#')
        {
            idKmer <<= one;
            idKmer |= ( (kmer & (one << (n - i - 1)) ) >> ( (n - i - 1) ) );
        }
        else
        {
            hashed <<= one;
            hashed |=  ( (kmer & (one << (n - i - 1)) ) >> ( (n - i - 1) ) ) ;
        }
    }
    
    node.rawKmer =  idKmer ;
    node.hashedKmer =  hashed ;
}


HashedNode YRJUnieative::convetToHashed(LONG kmer)
{
    HashedNode node;
    node.index = this->getIndex();
    extractKeys(node, kmer);
    return node;
}









void YRJUnieative::fillTheHashedNodesVector()
{
    //to constrcut the data
    this->fillTheKmersVector();
    
    this->hashedKmersSize = this->getNumOfKmers();
    this->hashedKmers.resize(this->hashedKmersSize);
    
    for(LONGS i = 0 , n = hashedKmers.size() ; i < n ; ++ i)
    {
        this->hashedKmers[i] = this->convetToHashed(this->kmersVector[i]);
    }
    
    this->clearTheCompleteKmers(); // deleting the old object data
}




