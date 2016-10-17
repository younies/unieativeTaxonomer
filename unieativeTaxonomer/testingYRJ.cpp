//
//  testingYRJ.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/15/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "testingYRJ.hpp"



void Tester::testYRJvector(){
    ifstream simBA5Stream(this->path_to_simBA5);
    ifstream headerSimBA5Stream(this->path_to_simBA5_headers);
    
    
    
    string header;
    string DNA;
    
    this->hashes.resize(6);
    this->hashes[0] = new Hash(pattern1 , path_to_the_hashed_databases);
    this->hashes[1] = new Hash(pattern11 , path_to_the_hashed_databases);
    this->hashes[2] = new Hash(pattern2 , path_to_the_hashed_databases);
    this->hashes[3] = new Hash(pattern3 , path_to_the_hashed_databases);
    this->hashes[4] = new Hash(pattern4 , path_to_the_hashed_databases);
    this->hashes[5] = new Hash(pattern5 , path_to_the_hashed_databases);
    
    
    
    
    ofstream * result = new ofstream(this->result);
    
    
    while (getline(simBA5Stream, DNA )) {
        getline(simBA5Stream, DNA );
        YRJObject yrj(DNA);
        this->test1YRJ(result, yrj);
    
    }
    
}




void Tester::test1YRJ(ofstream * result , YRJObject & yrj)
{
    
    for (int i = 1 ; i < 7 ; ++i)
    {
        auto hits = this->getHitsWithDifference(&yrj, 4, i);
        
        *result << hits.size() << "\t" ;
    }
    
    *result << endl;
    
}



set<short> Tester::getHitsWithDifference(YRJObject * yrj , int diff , int hashNumber){
    
    set<short> hits;
    
    for(auto kmer: yrj->kmersVector)
    {
        for(int i = 0 , n = (int)this->hashes.size() ; i < n && i < hashNumber ; ++i)
        {
            auto hash = this->hashes[i];
            auto hit_diff_pairs = this->getNumerOfDifferences(hash, kmer);
            
            for(auto hit : hit_diff_pairs)
                if(hit.second <= diff)
                    hits.insert(hit.first);
                
        }
    }
    
    return hits;
}













vector< pair<short, short> >  Tester::getNumerOfDifferences(Hash * hash , LONG kmer){
    vector< pair<short, short> >ret;
    
    //getting the streams
    ifstream &index = *(hash->getIndexStream());
    ifstream &data = *(hash->getDataStream());
    
    //getting the hashed Kmer
    auto convertedKmer = hash->getHashedNode(kmer);
    INT kmerIndex = hash->getINTfromPair(convertedKmer.rawKmer);
    index.seekg(kmerIndex  * sizeof(HashIndex));
    
    HashIndex hashIndex;
    index.read((char *) &hashIndex, sizeof(HashIndex));
    
    
    
    if(hashIndex.size == 0){
        cout << "NOT found for kmer :" << kmer << endl;
        return ret;
    }
    
    
    auto dataIndex = hashIndex.getIndex();
    
    data.seekg(dataIndex);
    
    vector<HashData> hashedData(hashIndex.size);
    
    for(int i = 0 , n = hashIndex.size ; i < n ; ++i)
    {
        data.read((char*) &hashedData[i], sizeof(HashData));
    }
    
    
    //add first kmer
    ret.emplace_back(make_pair(hashedData[0].index, this->numOfDifferencesBetweenKmers(hashedData[0].hashedKmer, convertedKmer.hashedKmer)));
    
    for(auto hashedBlock : hashedData){
        short diff = this->numOfDifferencesBetweenKmers(hashedBlock.hashedKmer, convertedKmer.hashedKmer);
        
        if(ret.back().first == hashedBlock.index)
            ret.back().second = min(ret.back().second , diff);
        else
            ret.emplace_back(make_pair(hashedBlock.index, diff));
    }
    
    
    return ret;
    
    
    
    
}






short Tester::numOfDifferencesBetweenKmers(pair<short, short> hashedKmer1 , pair<short, short> hashedKmer2)
{
    
    short differences = 0;
    
    
    while (hashedKmer1.first != hashedKmer2.first)
    {
        if(hashedKmer1.first % 4 != hashedKmer2.first % 4)
            ++differences;
        
        hashedKmer1.first /= 4;
        hashedKmer2.first /= 4;
    }
    
    while (hashedKmer1.second != hashedKmer2.second)
    {
        if(hashedKmer1.second % 4 != hashedKmer2.second % 4)
            ++differences;
        hashedKmer1.second /= 4;
        hashedKmer2.second /= 4;
    }
    
    return differences;
}


