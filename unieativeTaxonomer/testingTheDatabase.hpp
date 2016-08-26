//
//  testingTheDatabase.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/26/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef testingTheDatabase_hpp
#define testingTheDatabase_hpp


#include "headers.h"
#include "YRJObject.hpp"
#include "helpers.hpp"
#include "CoreTaxonomer.hpp"
#include "Tree.hpp"

void databaseTesting( string path_to_the_tree , string path_to_the_hashed_databases , string pattern , string path_to_the_yrj_databases , string path_to_million_random);

bool testKmer( CoreTaxonomer * core ,  Hash * hash, LONG kmer);

#endif /* testingTheDatabase_hpp */
