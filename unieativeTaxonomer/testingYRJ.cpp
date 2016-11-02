//
//  testingYRJ.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/15/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "testingYRJ.hpp"



void Tester::testYRJvector(int deep){
    ifstream simBA5Stream(this->pathMiSeq);
    ifstream headerSimBA5Stream(this->path_to_simBA5_headers);
    
    
    
    ifstream hiseqNew(this->hiSeqModified);
    
    
    
    this->pruinedTree = new Tree(this->path_to_pruined_tree);
    
    
    this->bigTree = new BigTree(path_to_the_names_dmp_file, path_to_the_nodes_dmp_file);
    
    string header;
    string DNA;
    
    this->hashes.resize(6);
    this->hashes[0] = new Hash(pattern1 , path_to_the_hashed_databases);
    this->hashes[1] = new Hash(pattern11 , path_to_the_hashed_databases);
    this->hashes[2] = new Hash(pattern2 , path_to_the_hashed_databases);
    this->hashes[3] = new Hash(pattern3 , path_to_the_hashed_databases);
    this->hashes[4] = new Hash(pattern4 , path_to_the_hashed_databases);
    this->hashes[5] = new Hash(pattern5 , path_to_the_hashed_databases);
    
    
    
    string s;
    s.push_back( (char)(deep + '0'));
    
    ofstream * result = new ofstream(this->result + s + ".out" );
    
    
    this->countHashEffect.resize(11, 0);
    this->globalCounter.resize(11,0);
    this->globalCounterEffect.resize(11, 0);
    
    
    /*
    vector<vector<string>> genomes;
    
    
    int i = -1;
    while (getline(simBA5Stream, DNA))
    {
        
        if(DNA[0] == '>'){
            genomes.emplace_back(vector<string>());
            ++i;
        }
        else
            genomes[i].emplace_back(DNA);
        
    }
    
    
    for ( auto genome : genomes) {
        
        YRJObject yrj(genome);
        if(yrj.getNumOfKmers() == 0){
            cout << "ooooh\n";
            continue;
            
        }
        
        this->testTheDifferences(result, yrj , 10);
    
    }
     
    */
    
    
    this->finalResult[this->notConsidered] = 0;
    this->finalResult[this->notNeeded] = 0;
    
    this->pruinedTree->bigTree = this->bigTree;
    
    while (getline(hiseqNew, DNA))
    {
        long uid = stol(DNA);
        getline(hiseqNew, DNA);
        YRJObject yrj(DNA);
        yrj.uid = uid;
        
        if (bigTree->getGenusUID(uid) == -1)
        {
            cout << "not  a good level " << uid << endl;
            continue;
        }
     //   cout << DNA << endl;
        this->testingGenomeLevel(&yrj, deep);
        //testingGenomeLevelWithNewMethodology(&yrj, deep);
    }
    
    
    calculate_accurcy(this->result + "_summery_" + s + ".out");
    
    
    string ranks[] = {"root" , "phylum" , "class" , "order" , "family" , "genus" , "species" , "subspecies" , "no" , notConsidered , notNeeded};
    
    for(auto rank : ranks)
    {
        *result << rank  << "\t" << getElementInTheResult(rank) << endl;
    }
    
    result->close();
    
    
    delete bigTree;
    delete pruinedTree;
    
    
}








void Tester::test1YRJ(ofstream * result , YRJObject & yrj)
{
    
    for (int i = 1 ; i < 7 ; ++i)
    {
        auto hits = this->getHitsWithDifference(&yrj, 4, i);
        
        *result << hits.size() << "\t" ;
        this->countHashEffect[i-1] += hits.size();
    }
    
    *result << endl;
    
    *result << endl;
    
}



void Tester::testTheDifferences(ofstream * result , YRJObject & yrj, int maxDiff)
{
    
    
    for(int diff = 0 ; diff <= maxDiff ; ++diff){
        
        auto hits = getHitsWithDifferenceButFullHahes(&yrj , diff);
        
        *result << hits.size() << '\t' ;
        globalCounter[diff] += hits.size();
        
        if(hits.size() == 0)
            countHashEffect[diff]++;
        else
            this->globalCounterEffect[diff]++;
    }
    
    *result << endl ;
    
    for(auto co : globalCounter)
        *result << co << '\t' ;
    *result << endl;
    
    for(auto co : countHashEffect)
        *result << co << '\t' ;
    *result << endl;
    
    for(auto co : globalCounterEffect)
        *result << co << '\t' ;
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





set<short> Tester::getHitsWithDifferenceButFullHahes(YRJObject * yrj , int diff ){
    
    set<short> hits;
    
    for(auto kmer: yrj->kmersVector)
    {
        for(auto hash : hashes)
        {
            auto tempHits = getNumerOfDifferences(hash, kmer);
            
            for( auto hit : tempHits)
            {
                if(hit.second <= diff)
                    hits.insert(hit.first);
            }
            
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
        //cout << "NOT found for kmer :" << kmer << endl;
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





vector< short > Tester::hits_kmer_with_differences( LONG kmer , int diff)
{
    set<short> ret_set;
    
    for (auto hash : this->hashes)
    {
        //extract the hits from each hash
        auto temp_hits = this->getNumerOfDifferences(hash, kmer);
        
        for (auto hit :  temp_hits)
        {
            if(hit.second <= diff)
                ret_set.insert(hit.first);
        }
    }
    
    vector<short> ret(ret_set.size()) ;
    
    int i = 0;
    
    for(auto retElement : ret_set)
    {
        ret[i++] = retElement;
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


