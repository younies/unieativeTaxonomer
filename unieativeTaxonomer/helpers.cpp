//
//  helpers.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/12/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "helpers.hpp"



/* Function to reverse bits of num */

LONG reverseKmer(LONG kmer)
{
    LONG  NO_OF_BITS = sizeof(kmer) * 8;
    LONG reverse_num = 0, i, temp;
    
    for (i = 0; i < NO_OF_BITS; i++)
    {
        temp = (kmer & (1 << i));
        if(temp)
            reverse_num |= (1 << ((NO_OF_BITS - 1) - i));
    }
    
    return reverse_num;
}
