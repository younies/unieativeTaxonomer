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
    
    output_file.write( (char *) &this->coreHashNodesSize ,  sizeof(LONG));
    
    for (LONGS i = 0 ; i < this->coreHashNodesSize ; ++i)
    {
        output_file.write( (char*)&this->coreHashedNodes[i] ,  sizeof(LONG));
    }
    
    output_file.flush();
    output_file.close();
}