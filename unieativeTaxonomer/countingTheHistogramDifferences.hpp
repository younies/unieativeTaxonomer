//
//  countingTheHistogramDifferences.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 9/11/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef countingTheHistogramDifferences_hpp
#define countingTheHistogramDifferences_hpp

#include "headers.h"
#include "YRJObject.hpp"
#include "helpers.hpp"
#include "CoreTaxonomer.hpp"
#include "Tree.hpp"

void calulateTheDifferences( string path_to_the_tree , string path_to_the_hashed_databases , string pattern , string path_to_the_yrj_databases , string path_to_million_random , string path_to_the_output_file);

vector<pair<short, int> > getNumberOfDifference( ifstream * indexStream , ifstream * dataStream,  Hash * hash, LONG kmer);

int numOfDifferencesBetweenKmers(pair<short, short> hashedKmer1 , pair<short, short> hashedKmer2);

#endif /* countingTheHistogramDifferences_hpp */
