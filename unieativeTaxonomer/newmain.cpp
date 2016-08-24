//
//  main.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "headers.h"
#include "YRJObject.hpp"
#include "helpers.hpp"
#include "CoreTaxonomer.hpp"
#include "Configurations.h"
#include "Tree.hpp"

string pattern = "##-#--###---#-#-#-#-#--#-##--##";
string pathe   = "/export1/project/hondius/newKrakenResearch/databases/kmerDatabase_new_31/all.yrj";
string million = "/export1/project/hondius/newKrakenResearch/generate_unique_random/MillionRandomIndices.txt";
string tenMillion = "/export1/project/hondius/newKrakenResearch/generate_unique_random/TenMillionRandom.txt";



string writeIN ="/export1/project/hondius/newKrakenResearch/generate_unique_random/";


vector<short> getIndiciesesAtDifferences(vector<pair< short , short> > & shortWithDifferences ,int maxDifferences);

bitset<64> updateHashValue(string hash);


HashedNode getTheHashedKmer(LONG kmer , bitset<64> hash);


LONG getTheKmer( HashedNode hashedNode , bitset<64> hash);





















































































bitset<64> updateHashValue(string hash)
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
    
    return newBitHash;
}














HashedNode getTheHashedKmer(LONG kmer , bitset<64> hash)
{
    HashedNode ret;
    
    bitset<sizeof(INT) * 8> first(0) , second(0);
    
    bitset<sizeof(kmer) * 8> kmerBits(kmer);
    
    int posFirst = 0 , posSecond = 0;
    for(LONGS i = 0 , n = hash.size() ; i <  n ; ++i )
    {
        if(hash[i])
        {
            first[posFirst++] = kmerBits[i];
        }
        else
        {
            second[posSecond++] = kmerBits[i];
        }
    }
    
    ret.index = -1;
    INT rawKmerINT = (INT)first.to_ulong();
    INT hahsedINT = (INT)second.to_ulong();
    ret.rawKmer.first = rawKmerINT >> (sizeof(SHORT) * 8);
    ret.rawKmer.second = (rawKmerINT << (sizeof(SHORT) * 8 ) ) >> (sizeof(SHORT) * 8);
    
    ret.hashedKmer.first = hahsedINT >> (sizeof(SHORT) * 8);
    ret.hashedKmer.second = (hahsedINT << (sizeof(SHORT) * 8 ) ) >> (sizeof(SHORT) * 8);
    
    
    return ret;
}



LONG getTheKmer( HashedNode hashedNode , bitset<64> hash)
{
    LONG ret;
    
    INT rawKmer = hashedNode.rawKmer.first;
    rawKmer <<= 16;
    rawKmer |= hashedNode.rawKmer.second;
    
    
    INT hashedKmer = hashedNode.hashedKmer.first;
    hashedKmer <<= 16;
    hashedKmer |= hashedNode.hashedKmer.second;

    bitset<32> bitRaw(rawKmer);
    bitset<32> bitHashed(hashedKmer);
    bitset<64> bitRet;
    
    
    int rawPos = 0 , hashPos = 0;
    
    for(int i = 0 ; i < 64 ; ++i)
    {
        if(hash[i])
            bitRet[i] = bitRaw[rawPos++];
        else
            bitRet[i] = bitHashed[hashPos++];
    }
    
    
    return bitRet.to_ulong();
}










int main(int argc, const char * argv[])
{
    int t = 1000;
    
    bitset<64> hash = updateHashValue(pattern);
    
    while(t--)
    {
        LONG rand = random();
        cout << rand << endl;
        if( rand == getTheKmer(getTheHashedKmer( rand,hash), hash) )
            cout << "YES\n";
        else
            cout << "NO\n";
    }
    
    return 0;
}





























































































































































































































































































































vector<short> getIndiciesesAtDifferences(vector<pair< short , short> > & shortWithDifferences ,int maxDifferences)
{
    vector<short> ret;
    
    for (LONGS i = 0 , n = shortWithDifferences.size() ; i < n ; ++i)
    {
        if(shortWithDifferences[i].second <= maxDifferences)
            ret.push_back(shortWithDifferences[i].first);
    }
    
    return ret;
}


