//
//  testingYRJ.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/15/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef testingYRJ_hpp
#define testingYRJ_hpp

#include "headers.h"
#include "YRJObject.hpp"
#include "BigTree.hpp"
#include "Tree.hpp"
#include "Hash.hpp"
#include "helpers.hpp"
#include "CoreTaxonomer.hpp"

class Tester {

    //constant we need in our program
    const string path_to_simBA5 = "/export1/project/hondius/testingUnieative/accuracy/simBA5_accuracy.fa";
    const string path_to_simBA5_headers = "/export1/project/hondius/testingUnieative/accuracy/mappingGI.txt";
    const string path_to_the_hashed_databases = "/export1/project/hondius/unieative/databases/";

    
    vector<YRJObject* > yrjVector;
    
    BigTree * bigTree;
    Tree * pruinedTree;
    
    vector<Hash *> hashes;
    

    
    
public:
    
    
    
    
    void testYRJvector();
    
    vector< pair<short, short> > getNumerOfDifferences(Hash * hash , LONG kmer);

    //return the difference between two hashed kmers
    short numOfDifferencesBetweenKmers(pair<short, short> hashedKmer1 , pair<short, short> hashedKmer2);

};


#endif /* testingYRJ_hpp */
