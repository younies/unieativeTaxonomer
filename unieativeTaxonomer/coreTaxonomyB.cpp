//
//  coreTaxonomyB.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/18/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "CoreTaxonomer.hpp"

void CoreTaxonomer::writeTheCoreData(string path)
{
    ofstream output_file(path , ios::binary);
    
    LONG hashLong = this->hash.to_ulong();
    output_file.write( (char *) &this->kmerLength ,  sizeof(LONG));
    output_file.write( (char *) &hashLong ,  sizeof(LONG));
    output_file.write( (char *) &this->coreHashNodesSize ,  sizeof(LONG));
    
    for (LONGS i = 0 ; i < this->coreHashNodesSize ; ++i)
    {
        output_file.write( (char*)&this->coreHashedNodes[i] ,  sizeof(LONG));
    }
    
    output_file.flush();
    output_file.close();
}


CoreTaxonomer::CoreTaxonomer(string hash , string path)
{
    //setting up the hash
    updateHashValue(hash);
    
    //
    
    ifstream fileStream(path);
    
    LONGS theSize;
    
    this->coreHashNodesSize = theSize;
    
    fileStream.read( (char *)&theSize  , sizeof(LONG));
    
    this->coreHashedNodes.resize(theSize);
    
    fileStream.read((char *) &this->coreHashedNodes[0], sizeof(HashedNode) * theSize);

    fileStream.close();
    
    
}





pair<LONGS, LONGS>  CoreTaxonomer::getThePlaceOfKmer(pair<SHORT, SHORT> rawKmer)
{
    pair<LONGS, LONGS> ret( -1 , -1 );
    
    LONGS start = 0 , end = this->coreHashNodesSize - 1;
    
    while (end >= start)
    {
        LONGS mid = (start + end) / 2;
        
        if(this->coreHashedNodes[mid].rawKmer == rawKmer)
        {
            ret.first = mid;
            end = mid - 1;
        }
        else if (this->coreHashedNodes[mid].rawKmer > rawKmer)
            end = mid - 1;
        else
            start = mid + 1;
    }
    
    
    start = 0 ;
    end = this->coreHashNodesSize - 1;
    
    while (end >= start)
    {
        LONGS mid = (start + end) / 2;
        if(this->coreHashedNodes[mid].rawKmer == rawKmer)
        {
            ret.second = mid;
            start = mid + 1;
        }
        else if (this->coreHashedNodes[mid].rawKmer > rawKmer)
            end = mid - 1;
        else
            start = mid + 1;
    }
    
    
    return ret;
}









/*

pair<LONGS, LONGS>  CoreTaxonomer::getThePlaceOfKmer(pair<SHORT, SHORT> rawKmer)
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
    
    if(start > end) return make_pair(-1, -1);
    
    LONGS newEndStart =  mid ,  newEndEnd  = mid;
    
    //for setting the start
    if(start >= 0 && start < this->coreHashNodesSize)
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
    if(end >= 0 && end < this->coreHashNodesSize)
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


*/




