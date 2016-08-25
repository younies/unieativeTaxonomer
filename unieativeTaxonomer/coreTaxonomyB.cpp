//
//  coreTaxonomyB.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/18/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "CoreTaxonomer.hpp"

void CoreTaxonomer::writeTheCoreData()
{
    
    ofstream index_output(this->theHash->getTheIndexPath() , ios::binary);
    ofstream data_output (this->theHash->getTheDataPath() , ios::binary);
    
    
    //fill the empty file
    HashIndex emptyIndex;
    for (LONGS i = 0 ; i < UINT_MAX; ++i)
    {
        index_output.write((char *) &emptyIndex , sizeof(emptyIndex));
    }
    
    index_output.flush();
    index_output.seekp(0); /// return to the start
    
    //initialize index file
    
    
    cout << "starting to write all the data ...." << endl;
    
    LONGS start = 0  , end = 0 , counter = 1;
    pair<SHORT, SHORT> rawIdentifier = this->coreHashedNodes[0].rawKmer;
    
    for (LONGS i = 1 ; i < this->coreHashNodesSize; ++i)
    {
        if(   this->coreHashedNodes[i].rawKmer == this->coreHashedNodes[i-1].rawKmer)
        {
            counter++;
            end++;
        }
        else
        {
            // now we need to write the data
            HashIndex hashIndex;
            hashIndex.size = counter;
            hashIndex.SetIndex(data_output.tellp());
            
            LONG index_in_the_index_file = rawIdentifier.first;
            index_in_the_index_file <<= 16;
            index_in_the_index_file |= rawIdentifier.second;
            
            index_output.seekp(index_in_the_index_file * sizeof(HashIndex));
            
            
            index_output.write((char *) &hashIndex , sizeof(HashIndex));
            
            
            for(LONGS j = start ; j <= end ; ++j)
            {
                HashData temp;
                temp.index = this->coreHashedNodes[j].index;
                temp.hashedKmer = this->coreHashedNodes[j].hashedKmer;
                data_output.write((char *) &temp, sizeof(temp));
            }
            
            // now we finished writing
            
            rawIdentifier = this->coreHashedNodes[i].rawKmer;
            counter = 1;
            start = i ; end = i;
            
        }
    }
    
    
    // now we need to write the data
    HashIndex hashIndex;
    hashIndex.size = counter;
    hashIndex.SetIndex(data_output.tellp());
    
    LONG index_in_the_index_file = rawIdentifier.first;
    index_in_the_index_file <<= 16;
    index_in_the_index_file |= rawIdentifier.second;
    
    index_output.seekp(index_in_the_index_file * sizeof(HashIndex));
    
    
    index_output.write((char *) &hashIndex , sizeof(HashIndex));
    
    
    for(LONGS j = start ; j <= end ; ++j)
    {
        HashData temp;
        temp.index = this->coreHashedNodes[j].index;
        temp.hashedKmer = this->coreHashedNodes[j].hashedKmer;
        data_output.write((char *) &temp, sizeof(temp));
    }
    
    
    
    

    
    data_output.flush();
    index_output.flush();
    
}


CoreTaxonomer::CoreTaxonomer(string hash , string path)
{
    //setting up the hash
    
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




