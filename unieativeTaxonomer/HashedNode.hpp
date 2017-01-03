//
//  HashedNode.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef HashedNode_hpp
#define HashedNode_hpp

#include "headers.h"
//#include "Hash.hpp"



/**
 
 This .hpp file include all the abstract for three types:
    1- HashedNode --> correspond to  a kmer after applying the hash
    2- HashData --> correspond to a kmer DataPart in the CoreDatabase
        -   SHortUID + Remaining Part
    3- HashIndex --> corresponding to a kmer Index part in the CoreDatabase.
        -   Index In The Data File + Size of Remaining Kmers
 
 */


//abstracting the new hashed kmer
struct HashedNode{
    short index;
    pair<SHORT, SHORT> rawKmer; //the non hashed part in the kmer
    pair<SHORT, SHORT> hashedKmer ; // the part that is hashed in the kmer
    
    HashedNode(  ){
        index = 0;
        rawKmer.first = 0;
        rawKmer.second = 0;
        hashedKmer.first = 0;
        hashedKmer.second = 0;
    }
    
};

bool  hashedNodeCompare( const HashedNode &lhs, const HashedNode &rhs);
//end of the hasehd kmer implementation


struct HashData
{
    short index;
    pair<SHORT, SHORT> hashedKmer ; // the part that is hashed in the kmer
    
    HashData(){
        index=0;
        hashedKmer.first = 0;
        hashedKmer.second = 0;
    }
    
};



struct HashIndex {
    pair<INT, INT> index;
    INT  size;
    
    HashIndex(){
        size = 0;
        index.first = 0;
        index.second = 0;
    }
    
    void SetIndex(LONG indx)
    {
        index.first = indx >> 32;
        index.second = (INT)indx ;
    }
    
    LONG getIndex()
    {
        LONG ret = 0;
        ret = index.first;
        ret <<= 32;
        ret |= index.second;
        
        return ret;
    }
    
};

#endif /* HashedNode_hpp */
