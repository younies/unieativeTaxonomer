//
//  YRJshortread.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "YRJObject.hpp"


YRJObject::YRJObject(vector<string> &shortReads  )
{
    this->kmerLength = this->kmerStandardLength;
    
    
    vector<string> reads(1, "");
    
    
    for(auto shortRead : shortReads)
    {
        for(auto c : shortRead)
        {
            if(c == 'N')
                reads.emplace_back("");
            else
                reads.back().push_back(c);
        }
    }
    
    this->fillKmersFromShortRead(reads);
    
    
    this->numOfKmers = this->kmersVector.size();
    
    
    
}




int YRJObject::getCorrespondingCode(char c)
{
    c = toupper(c);
    switch (c) {
        case 'A':
            return 0;
            break;
        case 'C':
            return 1;
            break;
        case 'G':
            return 2;
            break;
        case 'T':
            return 3;
            break;
            
        default:
            cout << "There are error\n";
            cout << c << endl;
            break;
    }
    
    return '\0';
}




LONG YRJObject::getLeastCanonicalKmer(const LONG kmer)
{
    //first inverse the kmer
    LONG invertKmer = kmer ^ this->maxLONG;
    invertKmer = invertKmer &  this->maxOnes31;
    
    
    //then revers it
    LONG canonicalKmer = 0;
    for (int i = 0 ;  i < this->kmerStandardLength ; ++i)
    {
        int letter = invertKmer % 4;
        canonicalKmer <<= 2;
        canonicalKmer |= letter;
        invertKmer >>= 2;
    }
    
    
    if(canonicalKmer < kmer)
        return canonicalKmer;
    return kmer;
    
}



void YRJObject::fillKmersFromShortRead(vector<string> & shortReads)
{
        
    unordered_set<LONG>  hashKmers;
    
    
    this->fillTheHashKmer(shortReads, hashKmers);
    
    //fill and sort the kmers' vector
        for(auto kmer : hashKmers)
            this->kmersVector.emplace_back(kmer);
        sort(this->kmersVector.begin(), this->kmersVector.end());
        this->numOfKmers = this->kmersVector.size();
    

}



void YRJObject::fillTheHashKmer(vector<string> & shortReads ,  unordered_set<LONG> & hashKmers)
{
    for(auto shortRead : shortReads)
    {
        if(shortRead.size() < this->kmerStandardLength)
            continue;
        LONG tempKmer = 0 ;
        //get the fisrt 30 characters in the short read
        for(int i = 0 ; i < this->kmerStandardLength - 1 ; ++i)
        {
            tempKmer <<= 2;
            int letter = this->getCorrespondingCode(shortRead[i]);
            tempKmer |= letter;
        }
        
        
        //go through all of the kmers
        for(int i = 30 , n = (int)shortRead.size() ; i < n ; ++i)
        {
            tempKmer <<= 2;
            tempKmer &= this->maxOnes31;
            int letter = this->getCorrespondingCode(shortRead[i]);
            tempKmer |= letter;
            
            hashKmers.insert(this->getLeastCanonicalKmer(tempKmer));
        }
    }
    
}



YRJObject::YRJObject(string &shortRead )
{
    
    
    vector<string> reads(1, "");
    
    
    
        for(auto c : shortRead)
        {
            if(c == 'N')
                reads.emplace_back("");
            else
                reads.back().push_back(c);
        }
        
    this->fillKmersFromShortRead(reads);
    
    
    this->numOfKmers = this->kmersVector.size();
    
}




