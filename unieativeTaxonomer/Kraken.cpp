//
//  Kraken.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/3/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "Kraken.hpp"



SHORT Unieative::getLCA(LONG kmer , int differences)
{
    
    //connectors to the datastreams and index streams
    vector<ifstream *>  indexStreams(this->Hashes.size());
    vector<ifstream *>  dataStreams(this->Hashes.size());
    
    
    
    
    for(auto hash: this->Hashes)
    {
        auto tempHits = getNumberOfDifference(hash->getTheIndexPath(), hash->getTheDataPath(), hash, kmer);
    }
}
