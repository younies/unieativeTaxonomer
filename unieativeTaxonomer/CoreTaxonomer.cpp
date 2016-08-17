//
//  CorseTaxonomer.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "CoreTaxonomer.hpp"


CoreTaxonomer::CoreTaxonomer(vector<YRJObject *> yrjVector , string hash)
{
    //setting up the hash
    updateHashValue(hash);
    
    this->yrjVector = yrjVector;
    //to calculate the size of the whole database hashed kmers
    this->coreHashNodesSize = 0;
    for(YRJObject* node : this->yrjVector)
        this->coreHashNodesSize += node->getNumOfKmers();
    
    //to build the hashed database
    this->coreHashedNodes.resize(this->coreHashNodesSize);
    
    //compy files in the coreHahsdNodes
    this->fillAllTheCoreData();
    
    
    //sort all the core data
    std::sort(this->coreHashedNodes.begin(), this->coreHashedNodes.end() , hashedNodeCompare);
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
   
    for ( YRJObject* yrj: this->yrjVector)
    {
        yrj->fillTheKmersVector();
        for (LONG kmer : yrj->kmersVector)
        {
            pair<INT, INT> tempHahsed =  getTheHashedKmer( kmer);
            this->coreHashedNodes[ this->startIndex].index = yrj->getIndex();
            this->coreHashedNodes[this->startIndex ].rawKmer = tempHahsed.first;
            this->coreHashedNodes[this->startIndex ].hashedKmer = tempHahsed.second;
            ++this->startIndex;
        }
    }
    
    if(this->startIndex == this->coreHashNodesSize)
    {
        cout << "perfect filling\n";
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
    while ( this->coreHashedNodes[start].rawKmer  != rawKmer)
    {
        ++start; // take a step to the goal
        mid = (start + newEndStart) / 2;
        
        if(this->coreHashedNodes[mid].rawKmer == rawKmer)
            newEndStart = mid;
        else
            start = mid + 1;
    }
    
    
    
    //for setting the end
    
    while (this->coreHashedNodes[end].rawKmer != rawKmer)
    {
        --end; //take one step to the goal and prevent the infinity loop!
        mid = (end + newEndEnd)/2;
        if(this->coreHashedNodes[mid].rawKmer == rawKmer)
            newEndEnd = mid;
        else
            end = mid  - 1;
    }
    
    return make_pair(start, end);
    
}




/**
 Implementing the function that convert the kmer from LONGS kmer to the corresponding rawKmer and hashedPart which is the hidden part
 */

pair<INT, INT> CoreTaxonomer::getTheHashedKmer(LONG kmer)
{
    pair<INT, INT> ret(0,0);
    
    bitset<sizeof(INT) * 8> first(0) , second(0);
    
    bitset<sizeof(kmer) * 8> kmerBits;
    
    int posFirst = 0 , posSecond = 0;
    for(LONGS i = 0 , n = this->hash.size() ; i <  n ; ++i )
    {
        if(this->hash[i])
        {
            first[posFirst++] = kmerBits[i];
        }
        else
        {
            second[posSecond++] = kmerBits[i];
        }
    }
    
    ret.first = (INT)first.to_ulong();
    ret.second = (INT)second.to_ulong();

    
    return ret;
}



void CoreTaxonomer::updateHashValue(string hash)
{
    //reverse the hash
    reverse(hash.begin() , hash.end());
    bitset<64> newBitHash;
    for(LONGS i = 0 , n = hash.size() ; i < n ; ++i)
    {
        LONGS ii = i * 2;
        if(hash[i] == '#' )
        {
            newBitHash[ii] = 1;
            newBitHash[ii+1] = 1;
        }
        else
        {
            newBitHash[ii] = 0;
            newBitHash[ii+1] = 0;
        }
    }
    
    this->hash = newBitHash;
}

/*
//reversing the bits of all the variable
LONG CoreTaxonomer::reverseKmer(LONG kmer)
{
    int size = sizeof(LONG) * 8;
    
    LONG ret = 0;
    int kmerLastBit = kmer%2;
    kmer /= 2;
    ret |= kmerLastBit;
    
    
    while (--size)
    {
        ret <<= 1;
        kmerLastBit = kmer%2;
        kmer /= 2;
        ret |= kmerLastBit;
    }
    
    
    return  ret;
}

*/

/**
 getting the short name from the the database
 */


vector<pair< short , short> > CoreTaxonomer::getShortNameFromKmer(LONG kmer)
{
    pair<INT, INT> hashedKmer = getTheHashedKmer(kmer);
    
    pair<LONGS, LONGS> startEnd = getThePlaceOfKmer(hashedKmer.first);
    

    vector< pair<short, short> > ret;
    ret.push_back(scanAtIndex(  startEnd.first , hashedKmer.second));

    LONGS curr = 0;
    for(LONGS i = startEnd.first + 1 ; i <= startEnd.second ; ++i)
    {
        pair<short, short> temp = scanAtIndex(i, hashedKmer.second);
        
        if( temp.first == ret[curr].first)
            ret[curr].second = min(ret[curr].second , temp.second );
        else
        {
            ++curr;
            ret.push_back(temp);
        }
        
    }
    
    return ret;
}



/**
 scanning the kmers

 */


pair<short, short>  CoreTaxonomer::scanAtIndex( LONGS index , INT hashed)
{
    pair<short, short> ret;
    ret.first = this->coreHashedNodes[index].index;
    ret.second = 0;
    
    INT comparedHash = this->coreHashedNodes[index].hashedKmer;
    
    int iterations = sizeof(INT) * 8;
    
    while (iterations --)
    {
        if(comparedHash%2 != hashed %2)
            ++ret.second;
        comparedHash /= 2;
        hashed /= 2;
    }
    
    return ret;
}






