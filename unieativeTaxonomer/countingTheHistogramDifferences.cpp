//
//  countingTheHistogramDifferences.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 9/11/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "countingTheHistogramDifferences.hpp"


//
//  testingTheDatabase.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/26/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "testingTheDatabase.hpp"


void calulateTheDifferences( string path_to_the_tree , string path_to_the_hashed_databases , string pattern , string path_to_the_yrj_databases , string path_to_million_random , string path_to_the_output_file)
{
    //random kmers
    YRJObject * randomMillionKmers = new YRJObject(path_to_million_random , -1);
    randomMillionKmers->fillTheKmersVector();
    
    //creating a tree object
    Tree * tree = new Tree(path_to_the_tree);
    
    //creating the hash
    Hash * hash = new Hash( pattern , path_to_the_hashed_databases );
    
    
    //now we want to have the connectors
    ifstream * indexStream = new ifstream(hash->getTheIndexPath());
    ifstream * dataStream  = new ifstream(hash->getTheDataPath());
    
    //now creat the output file
    ofstream * outputFileResults = new ofstream(path_to_the_output_file);
    
    
    vector<LONGS> globalCounter(7 , 0);
    
    for(LONGS i = 0 , n = randomMillionKmers->getNumOfKmers() ; i < n ; ++i)
    {
        vector<LONGS> localCounter(7,0);
        
        vector<pair<short, int> > tempResult = getNumberOfDifference(indexStream, dataStream, hash, randomMillionKmers->kmersVector[i]);
        
        for( pair<short, int> tempPair : tempResult)
        {
            if(tempPair.second > 5)
                localCounter[6]++;
            else
                localCounter[tempPair.second]++;
        }
        
    
        *outputFileResults << randomMillionKmers->kmersVector[i] << endl;
        for(int j = 0  , m = (int)localCounter.size() ; j < m ; ++j)
        {
            *outputFileResults << localCounter[j] << "\t" ;
            globalCounter[j] += localCounter[j];
            
        }
        *outputFileResults << endl;
        
        for(int j = 0  , m = (int)localCounter.size() ; j < m ; ++j)
            *outputFileResults << globalCounter[j] << "\t";
        *outputFileResults << endl;
        
        
    }
    
    outputFileResults->close();
    dataStream ->close();
    indexStream->close();
    
}






vector<pair<short, int> > getNumberOfDifference( ifstream * indexStream , ifstream * dataStream,  Hash * hash, LONG kmer)
{

    
    //ifstream indexStream(hash->getTheIndexPath());
    //ifstream dataStream (hash->getTheDataPath()) ;
    
    HashedNode kmerNode = hash->getHashedNode(kmer);
    
    INT index = hash->getINTfromPair(kmerNode.rawKmer);
    
    
    indexStream->seekg(index * sizeof(HashIndex));
    
    HashIndex tempHashIndex;
    
    indexStream->read((char *) &tempHashIndex, sizeof(HashIndex));
    
    
    LONG dataIndex = tempHashIndex.getIndex();
    
    dataStream->seekg(dataIndex );
    
    vector<HashData> hashDataVector(tempHashIndex.size);
    
    for(LONGS i = 0 ; i <tempHashIndex.size ; ++i)
    {
        dataStream->read((char *) &hashDataVector[i], sizeof(HashData));
    }
    
    vector<pair<short, int> > ret;
    
    for (int i = 0 , n = (int)hashDataVector.size(); i < n ; ++i)
    {
        if(ret.size() == 0)
            ret.push_back(make_pair(hashDataVector[i].index,  numOfDifferencesBetweenKmers(hashDataVector[i].hashedKmer, kmerNode.hashedKmer) ));
            
        else
        {
            int temDiff = numOfDifferencesBetweenKmers(hashDataVector[i].hashedKmer, kmerNode.hashedKmer);
            if(ret.back().first != hashDataVector[i].index )
            {
                ret.push_back(make_pair(hashDataVector[i].index, temDiff));
            }
            else
            {
                ret.back().second = min(temDiff , ret.back().second);
            }
        }
    }
    
    
    return ret;
}




int numOfDifferencesBetweenKmers(pair<short, short> hashedKmer1 , pair<short, short> hashedKmer2)
{
    
    int differences = 0;
    
    
    while (hashedKmer1.first != hashedKmer2.first)
    {
        differences++;
        hashedKmer1.first /= 4;
        hashedKmer2.first /= 4;
    }
    
    while (hashedKmer1.second != hashedKmer2.second)
    {
        differences++;
        hashedKmer1.second /= 4;
        hashedKmer2.second /= 4;
    }
    
    return differences;
}




