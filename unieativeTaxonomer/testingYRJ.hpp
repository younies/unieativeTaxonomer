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

    
    //conf main
    const string pattern1 =  "##-#--###---#-#-#-#-#--#-##--##";
    const string pattern11 = "#-#-##---#-#-#-#-#-#-##-#--##-#";
    const string pattern2 = "##-##-#-#---#-#-##--#--#-#-#-##";
    const string pattern3 = "#-#-----#-#--#--##-##-#-####-##";
    const string pattern4 = "#--###---#--##--#-#-###-#-##--#";
    const string pattern5 = "#-####-##----#-##---#-#-#--##-#";
    
    const string result = "/export1/project/hondius/testingUnieative/results/thesis_trial2.txt";
    vector<YRJObject* > yrjVector;
    
    BigTree * bigTree;
    Tree * pruinedTree;
    
    vector<Hash *> hashes;
    

    vector<LONGS> countHashEffect;
    
    
public:
    
    
    
    
    void testYRJvector();
    
    vector< pair<short, short> > getNumerOfDifferences(Hash * hash , LONG kmer);

    //return the difference between two hashed kmers
    short numOfDifferencesBetweenKmers(pair<short, short> hashedKmer1 , pair<short, short> hashedKmer2);

    void test1YRJ(ofstream * result , YRJObject & yrj);
    set<short> getHitsWithDifference(YRJObject * yrj , int diff , int hashNumber);
    
};


#endif /* testingYRJ_hpp */
