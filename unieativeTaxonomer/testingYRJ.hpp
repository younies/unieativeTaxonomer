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
    
    const string pathHiSeq = "/export1/project/hondius/testingUnieative/accuracy/HiSeq_accuracy.fa";
    const string pathMiSeq = "/export1/project/hondius/testingUnieative/accuracy/MiSeq_accuracy.fa";
    const string path_to_simBA5 = "/export1/project/hondius/testingUnieative/accuracy/simBA5_accuracy.fa";
    const string path_to_simBA5_headers = "/export1/project/hondius/testingUnieative/accuracy/mappingGI.txt";
    const string path_to_the_hashed_databases = "/export1/project/hondius/unieative/databases/";

    const string path_to_pruined_tree = "/export1/project/hondius/newKrakenResearch/finalNewTree.txt";

    //conf main
    const string pattern1 =  "##-#--###---#-#-#-#-#--#-##--##";
    const string pattern11 = "#-#-##---#-#-#-#-#-#-##-#--##-#";
    const string pattern2 = "##-##-#-#---#-#-##--#--#-#-#-##";
    const string pattern3 = "#-#-----#-#--#--##-##-#-####-##";
    const string pattern4 = "#--###---#--##--#-#-###-#-##--#";
    const string pattern5 = "#-####-##----#-##---#-#-#--##-#";
    const string path_to_the_names_dmp_file = "/export1/project/hondius/newKrakenResearch/databases/names.txt";
    const string path_to_the_nodes_dmp_file = "/export1/project/hondius/newKrakenResearch/databases/nodes.txt";

    const string result = "/export1/project/hondius/testingUnieative/results/MiSeqMaxDiff.txt";
    vector<YRJObject* > yrjVector;
    
    BigTree * bigTree;
    Tree * pruinedTree;
    
    unordered_map<string , long> finalResult;
    
    vector<Hash *> hashes;
    

    vector<LONGS> countHashEffect;
    vector<long> globalCounter;
    vector<long> globalCounterEffect;
    
public:
    
    
    
    
    void testYRJvector();
    
    vector< short > hits_kmer_with_differences( LONG kmer , int diff);// return all the hits

    
    
    vector< pair<short, short> > getNumerOfDifferences(Hash * hash , LONG kmer);

    //return the difference between two hashed kmers
    short numOfDifferencesBetweenKmers(pair<short, short> hashedKmer1 , pair<short, short> hashedKmer2);

    void test1YRJ(ofstream * result , YRJObject & yrj);
    set<short> getHitsWithDifference(YRJObject * yrj , int diff , int hashNumber);
    
    set<short> getHitsWithDifferenceButFullHahes(YRJObject * yrj , int diff );
    
    void testTheDifferences(ofstream * result , YRJObject & yrj , int maxDiff);
    
    void testingGenomeLevel( YRJObject * yrj  ,int differences);
};


#endif /* testingYRJ_hpp */
