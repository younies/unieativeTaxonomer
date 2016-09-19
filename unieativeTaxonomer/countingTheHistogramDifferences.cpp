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


void CalculateDifferencesDistributions(  string path_to_the_hashed_databases , vector<string> patterns , string path_to_the_yrj_databases , string path_to_million_random , string path_to_the_output_file)

{
    ofstream * writtingFile = new ofstream(path_to_the_output_file);
    
    int numberOfHashes = (int)patterns.size();
    
    //random kmers
    YRJObject * randomMillionKmers = new YRJObject(path_to_million_random , -1);
    randomMillionKmers->fillTheKmersVector();
    
    //create the hash vector
    vector<Hash *> hashes(patterns.size());
    for (int i = 0 , n = (int)hashes.size() ; i < n ; ++i)
    {
        hashes[i] = new Hash(patterns[i] , path_to_the_hashed_databases);
    }
    
    //connector to data files
    vector<ifstream *>  indexStreams(patterns.size());
    vector<ifstream *>  dataStreams(patterns.size());
    
    for (int i = 0 , n = (int)hashes.size(); i < n  ; ++i)
    {
        indexStreams[i] = new ifstream(hashes[i]->getTheIndexPath());
        dataStreams[i]  = new ifstream(hashes[i]->getTheDataPath());
    }
    
    //open the writing files
    vector<ofstream *> outputFilesResults(numberLimit);
    for (int i =  0  ;  i < numberLimit ; ++i)
    {
        outputFilesResults[i ] = new ofstream(path_to_the_output_file + "_" + to_string(i) + ".result");
    }
    
    vector< vector< LONGS> > globalCounter( numberOfHashes , vector<LONGS>( numberLimit , 0) );
    
    
    for(LONGS i = 0 , n = randomMillionKmers->getNumOfKmers() ; i < n ; ++i)
    {
         vector<vector< vector<SHORT> > > localCounter(numberOfHashes, vector<vector<SHORT> > (numberLimit));
        cout << "  hashes  " << hashes.size();
        cout << i << "  " << n << endl;
        for (LONGS local_i , local_n = hashes.size();  local_i < local_n; ++local_i )
        {
            vector<pair<short, int> > tempResult = getNumberOfDifference(indexStreams[local_i], dataStreams[local_i], hashes[local_i], randomMillionKmers->kmersVector[i]);

            for( auto tempPair : tempResult)
            {
                if(tempPair.second >= numberLimit - 1)
                    localCounter[local_i][numberLimit - 1].emplace_back(tempPair.first);
                else
                    localCounter[local_i][tempPair.second].emplace_back(tempPair.first);
            }

        }
        
        ///Merge Tomorrwo for a random kmer  [i]
        for(LONGS hash_i =  1  , hash_n = hashes.size() ; hash_i < hash_n ; ++hash_i )
            mergeTwoColumnsToColumnTwo(localCounter[hash_i - 1], localCounter[hash_i]);
        
        
        //preparing for writing
        vector<vector<LONGS>> tempLocal =  addLocalToGlobal(globalCounter , localCounter);
        
        *writtingFile << "Local " << i + 1 << endl;
        writeMatrixInFile( writtingFile , tempLocal);
        
        *writtingFile << "global " << i + 1 << endl;
        writeMatrixInFile( writtingFile , globalCounter);
        
        
        
        
        writtingFile->flush();
        
    }
 
    writtingFile->close();
}



void calulateTheDifferences( string path_to_the_tree , string path_to_the_hashed_databases , string pattern , string path_to_the_yrj_databases , string path_to_million_random , string path_to_the_output_file)
{
    //random kmers
    YRJObject * randomMillionKmers = new YRJObject(path_to_million_random , -1);
    randomMillionKmers->fillTheKmersVector();

    
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
    cout << "called " << endl;

    
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
    
    //testing
    //for(auto tempair: ret)
      //  cout << tempair.first << " " << tempair.second << endl;
    
    
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




vector<vector<short>> getDifferencesVector(vector<pair<short, int> >  differences)
{
    vector<vector<short>>  ret(numberLimit);
    
    for(pair<short, int> element: differences)
    {
        if(element.second < numberLimit - 1)
            ret[element.second].push_back(element.first);
        else
            ret[numberLimit - 1].push_back(element.first);
    }
    
    return ret;
}




void mergeTwoColumnsToColumnTwo(const vector< vector<SHORT> >&  vec1 , vector< vector<SHORT> >&  vec2)
{
    for (int i = 0 , n = (int) vec1.size() ; i < n  ; ++i)
    {
        vec2[i] = mergeTwoWithoutRepetition(vec1[i], vec2[i]);
    }
}


vector<SHORT> mergeTwoWithoutRepetition(const vector<SHORT> & vec1 , const vector<SHORT> & vec2)
{
    vector<SHORT> ret;
    LONGS retSize = 0;
    
    LONGS n1 = vec1.size()  , n2 = vec2.size();
    LONGS i1 = 0 , i2 = 0;
    while (i1 < n1 || i2 < n2)
    {
        short value;
        if( i1 < n1 && i2 < n2  )
        {
            if(vec1[i1] < vec2[i2])
                value = vec1[i1++];
            else if (vec1[i1] > vec2[i2])
                value = vec2[i2++];
            else//both of them equal
            {
                value = vec1[i1++];
                ++i2;
            }
        }
        else if(i1 < n1)
            value = vec1[i1++];
        else
            value = vec2[i2++];
        
        if(retSize == 0 || ret.back() != value )
        {
            ret.emplace_back(value);
            retSize++;
        }
    }
    
    return ret;
}
















