//
//  Hash.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/24/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//


/**
 - Hash class is reposible for all the operation including the hashing, from converting a kmer to two hash part and remaining part, writing a kmer, compare two remaining parts of a kmer and return the hashed kmer to the original one
 
 - So, the Hash object needs to take only a string hash as a parameter (string of 16_'#' and 15_'-'  ), and the Hash class are going to do the rest.

 
 */

#ifndef Hash_hpp
#define Hash_hpp

#include "headers.h"
#include "HashedNode.hpp"


class Hash
{
    int  bitSize = 64;
    const string INDEX_EXTENTION = ".hash_index";
    const string DATA_EXTENSTION = ".hash_data";
    string hash;
    string name;
    string path_to_index_file;
    string path_to_data_file;
    
    ifstream * indexStream;
    ifstream * dataStream;
    
    bitset<64> bitHash;
    
public:
    Hash( string hash , string path_to_database  );
    
    HashedNode getHashedNode( LONG kmer);
    
    
    pair<SHORT, SHORT> convertPartBitSetToPair(bitset<32> bitKmerPart);
    pair<SHORT, SHORT> convertINTPartToPair(INT kmerPart);
    
    INT getINTfromPair(pair<SHORT, SHORT> partKmer);
    
    string getTheIndexPath(){return this->path_to_index_file;}
    string getTheDataPath(){return this->path_to_data_file;}
    
    
    void openConnector();
    void closeConnectors();
    
    ifstream * getIndexStream();
    ifstream * getDataStream();
    
    
    
};


string hash_check(string hash);// checking the validity of the hash and return the extended hash

#endif /* Hash_hpp */
