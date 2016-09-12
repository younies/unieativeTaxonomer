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
    
    
    //now we want to test
    
    LONGS countTrue = 0 , countFalse = 0;
    
    for(LONGS i = 0 , n = randomMillionKmers->getNumOfKmers() ; i < n ; ++i)
    {
        if( testKmer(core, hash, randomMillionKmers->kmersVector[i]))
            countTrue++;
        else
            countFalse++;
        
        cout << "CountTrue = " << countTrue << "   countFalse =  " << countFalse << endl;
    }
    
    
}






bool testKmer( CoreTaxonomer * core ,  Hash * hash, LONG kmer)
{
    bool returnflag = true;
    ifstream indexStream(hash->getTheIndexPath());
    ifstream dataStream (hash->getTheDataPath()) ;
    
    HashedNode node = hash->getHashedNode(kmer);
    pair<LONGS, LONGS> startEnd = core->getThePlaceOfKmer(node.rawKmer);
    
    INT index = hash->getINTfromPair(node.rawKmer);
    
    
    indexStream.seekg(index * sizeof(HashIndex));
    
    HashIndex tempHashIndex;
    
    indexStream.read((char *) &tempHashIndex, sizeof(HashIndex));
    
    
    if( tempHashIndex.size != (startEnd.second - startEnd.first + 1)  )
    {
        cout << "the sizes does not matsh\n";
        cout << "for kmer " + to_string(kmer) << endl;
        cout << "the hash size " << tempHashIndex.size << endl;
        cout << " the index returned is " << startEnd.first << "  " << startEnd.second << endl;
        returnflag = false;
    }
    
    
    LONG dataIndex = tempHashIndex.getIndex();
    
    dataStream.seekg(dataIndex );
    
    
    for(LONGS i = startEnd.first ; i <= startEnd.second ; ++i)
    {
        HashedNode completedHashNode = core->getHashedNode(i);
        HashData tempHashData;
        dataStream.read((char *) &tempHashData, sizeof(HashData));
        
        
        
        
        
        cout << "the shorts does  matsh " << tempHashData.index << "  " <<completedHashNode.index << endl;
        cout << "the hashed kmers does equal to each others" << endl;
        cout << "hash one " <<tempHashData.hashedKmer.first << "  " << tempHashData.hashedKmer.second;
        cout << "  second  " << completedHashNode.hashedKmer.first << "  " << completedHashNode.hashedKmer.second << endl;
        
        
        
        
        if(tempHashData.index != completedHashNode.index)
        {
            cout << "the shorts does not matsh " << tempHashData.index << "  " <<completedHashNode.index << endl;
            cout  << " kmer  " << kmer << endl;
        }
        
        if(tempHashData.hashedKmer != completedHashNode.hashedKmer)
        {
            cout << "the hashed kmers does not equal to each others" << endl;
            cout << "hash one " <<tempHashData.hashedKmer.first << "  " << tempHashData.hashedKmer.second;
            cout << "  second  " << completedHashNode.hashedKmer.first << "  " << completedHashNode.hashedKmer.second << endl;
            
            returnflag = false;
        }
        
        
    }
    
    return true;
}









