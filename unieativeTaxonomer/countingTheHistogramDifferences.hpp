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


#define numberLimit 14


/*
 This class was for providing statistics that we used in our report such as
        - How  many hash that we need to find a specific amount of hitted kmers
        - How many differences in the in the hitted kmers
        ..etc.
 */

vector<vector<short>> getDifferencesVector(vector<pair<short, int> >  differences);



// calculating the percentage in 1000,000 kmers. How many of them hit with zero differences, 1 differences ... etc.
void calulateTheDifferences( string path_to_the_tree , string path_to_the_hashed_databases , string pattern , string path_to_the_yrj_databases , string path_to_million_random , string path_to_the_output_file);

//for a specific hash, and specific kmer ... it will return the number of differences that this kmer hit with.
vector<pair<short, int> > getNumberOfDifference( ifstream * indexStream , ifstream * dataStream,  Hash * hash, LONG kmer);

//return the differences between two hashed kmers
int numOfDifferencesBetweenKmers(pair<short, short> hashedKmer1 , pair<short, short> hashedKmer2);


// this function merge two (sorted) vectors into one merged vector without repition
vector<SHORT> mergeTwoWithoutRepetition(const vector<SHORT> & vec1 , const vector<SHORT> & vec2);


// this function merge two matrix and no repition on columns vectors.
void mergeTwoColumnsToColumnTwo(const vector< vector<SHORT> >&  vec1 , vector< vector<SHORT> >&  vec2);


// write in the output file (path) how many kmers hit with 0 difference, 1 difference and so on.
void CalculateDifferencesDistributions(  string path_to_the_hashed_databases , vector<string> patterns , string path_to_the_yrj_databases , string path_to_million_random , string path_to_the_output_file);


void writeMatrixInFile(ofstream * writtingFile , vector<vector<LONGS> > data);


vector<vector< LONGS> >  addLocalToGlobal( vector< vector< LONGS> > &globalCounter , const vector<vector<vector<SHORT>>> & localCounter);




#endif /* countingTheHistogramDifferences_hpp */
