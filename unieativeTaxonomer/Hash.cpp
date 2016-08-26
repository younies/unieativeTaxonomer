//
//  Hash.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/24/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "Hash.hpp"



Hash::Hash( string hash , string path_to_database  )
{
    this->hash = hash_check(hash);
    
    string converterHash = hash_check(hash);
    
    reverse(converterHash.begin(), converterHash.end());
    
    for(int i = 0  , n = (int) converterHash.size() ; i < n ; ++i)
    {
        if( converterHash[i] == '#')
            this->bitHash[i] = 1;
        else
            this->bitHash[i] = 0;
    }//the bithash are sitted
    
    
    //set the name
    this->name = to_string(bitHash.to_ulong());
    
    //set the pathes
    this->path_to_index_file = path_to_database + name + this->INDEX_EXTENTION;
    this->path_to_data_file  = path_to_database + name + this->DATA_EXTENSTION;
    
    
    
}







//implementing the get hash function

HashedNode  Hash::getHashedNode( LONG kmer)
{
    bitset< 64 > bitKmer(kmer);
    HashedNode ret ;
    
    int rawPos = 0 , hashedPos = 0;
    bitset<32> rawPart , hashedPart;
    for(int i = 0 , n = (int) this->bitHash.size() ; i < n ; ++i)
    {
        if(bitHash[i])
            rawPart[rawPos++] = bitKmer[i];
        else
            hashedPart[hashedPos++] = bitKmer[i];
    }
    
    ret.rawKmer = convertPartBitSetToPair(rawPart);
    ret.hashedKmer = convertPartBitSetToPair(hashedPart);
    
    return ret;
}



pair<SHORT, SHORT> Hash::convertPartBitSetToPair(bitset<32> bitKmerPart)
{
    INT kmerPart = (INT)bitKmerPart.to_ulong();
    return this->convertINTPartToPair(kmerPart);
}


pair<SHORT, SHORT> Hash::convertINTPartToPair(INT kmerPart)
{
    pair<SHORT, SHORT> ret;
    
    ret.first = (kmerPart >> 16) ;
    ret.second = ( kmerPart << 16) >> 16 ;
    
    return ret;
}


string hash_check(string hash)
{
    int chbang = 0 , dag = 0  , others = 0;
    
    vector<char> retHash;
    
    for(int i = 0 , n= (int)hash.size() ; i < n ; ++i)//to count the number of 3, - , and others
    {
        switch (hash[i])
        {
            case '#':
                chbang++;
                break;
            case '-':
                dag++;
                break;
                
            default:
                others++;
                break;
        }
        
        retHash.push_back(hash[i]);
        retHash.push_back(hash[i]);
    }
    
    
    if(chbang == 16 &&  dag == 15 && others == 0)
    {
        string ret(retHash.begin() , retHash.end());
        return  ret;
    }
    else
    {
        cout << " the hash is not satified the conditions\n";
        string ret(retHash.begin() , retHash.end());
        return  ret;
    }
    
}





INT Hash::getINTfromPair(pair<SHORT, SHORT> partKmer)
{
    INT ret = 0 ;
    
    ret = partKmer.first;
    ret <<= 16;
    ret |= partKmer.second;
    
    return ret;
    
}
