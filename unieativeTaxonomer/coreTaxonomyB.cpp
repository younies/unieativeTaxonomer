//
//  coreTaxonomyB.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/18/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "CoreTaxonomer.hpp"

void CoreTaxonomer::writeTheCoreData(string path)
{
    ofstream output_file(path , ios::binary);
    
    LONG hashLong = this->hash.to_ulong();
    output_file.write( (char *) &this->kmerLength ,  sizeof(LONG));
    output_file.write( (char *) &hashLong ,  sizeof(LONG));
    output_file.write( (char *) &this->coreHashNodesSize ,  sizeof(LONG));
    
    for (LONGS i = 0 ; i < this->coreHashNodesSize ; ++i)
    {
        output_file.write( (char*)&this->coreHashedNodes[i] ,  sizeof(LONG));
    }
    
    output_file.flush();
    output_file.close();
}


CoreTaxonomer::CoreTaxonomer(string hash , string path)
{
    //setting up the hash
    updateHashValue(hash);
    
    //
    
    ifstream fileStream(path);
    
    LONGS theSize;
    
    this->coreHashNodesSize = theSize;
    
    fileStream.read( (char *)&theSize  , sizeof(LONG));
    
    this->coreHashedNodes.resize(theSize);
    
    fileStream.read((char *) &this->coreHashedNodes[0], sizeof(HashedNode) * theSize);

    fileStream.close();
    
    
}
