//
//  Hash.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/24/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef Hash_hpp
#define Hash_hpp

#include "headers.h"
#include "HashedNode.hpp"


class Hash
{
    int  bitSize = 64;
    const string INDEX_EXTENTION = ".hash_index";
    const string DATA_EXTENSTION = ".hash_index";
    string hash;
    string name;
    string path_to_index_file;
    string path_to_data_file;
    
    bitset<64> bitHash;
    
public:
    Hash( string hash , string path_to_database  );
    
    HashedNode getHashedNode( LONG kmer);
    
    
    pair<SHORT, SHORT> convertPartBitSetToPair(bitset<32> bitKmerPart);
    pair<SHORT, SHORT> convertINTPartToPair(INT kmerPart);
    
    
    
    
};


string hash_check(string hash);// checking the validity of the hash and return the extended hash

#endif /* Hash_hpp */
