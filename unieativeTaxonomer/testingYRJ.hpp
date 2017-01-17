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

/*
 This calss has the most essential engines that use Kraken and other implementation to identify
 the corresponding UID for an input short read and also to identify the result.
 */

class Tester {

    //All the constants is pathes to required files in our program
    //that are going to be readed from the hard drive
    
    //constant we need in our program
    
    const string pathHiSeq = "/export1/project/hondius/testingUnieative/accuracy/HiSeq_accuracy.fa";
    const string pathMiSeq = "/export1/project/hondius/testingUnieative/accuracy/MiSeq_accuracy.fa";
    const string path_to_simBA5 = "/export1/project/hondius/testingUnieative/accuracy/simBA5_accuracy.fa";
    const string path_to_simBA5_headers = "/export1/project/hondius/testingUnieative/accuracy/mappingGI.txt";
    const string path_to_the_hashed_databases = "/export1/project/hondius/unieative/databases/";

    const string path_to_pruined_tree = "/export1/project/hondius/newKrakenResearch/finalNewTree.txt";
    const string hiSeqModified = "/export1/project/hondius/testingUnieative/accuracy/hiSeqModified.fa";
    const string mi_SeqModified ="/export1/project/hondius/testingUnieative/accuracy/mi2.txt";

    const string notConsidered = "not_considered";
    const string notNeeded = "not_needed";
    //conf main
    const string pattern1 =  "##-#--###---#-#-#-#-#--#-##--##";
    const string pattern11 = "#-#-##---#-#-#-#-#-#-##-#--##-#";
    const string pattern2 = "##-##-#-#---#-#-##--#--#-#-#-##";
    const string pattern3 = "#-#-----#-#--#--##-##-#-####-##";
    const string pattern4 = "#--###---#--##--#-#-###-#-##--#";
    const string pattern5 = "#-####-##----#-##---#-#-#--##-#";
    const string path_to_the_names_dmp_file = "/export1/project/hondius/newKrakenResearch/databases/names.txt";
    const string path_to_the_nodes_dmp_file = "/export1/project/hondius/newKrakenResearch/databases/nodes.txt";

    const string result = "/export1/project/hondius/testingUnieative/newNewTest/krakenPure/kmer_distribution_mi_mi_";
    
    const string path_to_Kraken_test = "/export1/project/hondius/testingUnieative/finalResults/ids_for_kraken_HiSeq.txt";
    
    vector<YRJObject* > yrjVector;
    
    BigTree * bigTree;
    Tree * pruinedTree;
    
    
    vector<unordered_map<string , long> > finalKmerResult;
    
    vector<pair<unordered_map<string , long> , unordered_map<string , long> >> max_min_final_results;
    
    unordered_map<string , long> finalResult;
    
    
    unordered_map<string , long> levelsValues;
    
    vector<Hash *> hashes;
    

    vector<LONGS> countHashEffect;
    vector<long> globalCounter;
    vector<long> globalCounterEffect;
    
    vector<vector<LONG> > kmers_block;
    
public:
    
    
    
    
    void testYRJvector(int deep);
    
    vector< short > hits_kmer_with_differences( LONG kmer , int diff);// return all the hits

    
    
    vector< pair<short, short> > getNumerOfDifferences(Hash * hash , LONG kmer);

    //return the difference between two hashed kmers
    short numOfDifferencesBetweenKmers(pair<short, short> hashedKmer1 , pair<short, short> hashedKmer2);

    void test1YRJ(ofstream * result , YRJObject & yrj);
    set<short> getHitsWithDifference(YRJObject * yrj , int diff , int hashNumber);
    
    set<short> getHitsWithDifferenceButFullHahes(YRJObject * yrj , int diff );
    
    void testTheDifferences(ofstream * result , YRJObject & yrj , int maxDiff);
    
    void testingGenomeLevel( YRJObject * yrj  ,int differences);
    
    // this function test if there is any hit at zero level?
    bool isKrakenCatch(YRJObject * yrj );
    // this function test if there is any hit at deep level?
    bool isKrakenCatch(YRJObject * yrj , int deep );
    
    
    map<short, int> getKrakenLCAs(YRJObject * yrj , int differences);
    
    map<LONGS, int>  getUnieativeHitsSpecies(YRJObject * yrj , int differences );

    // this method used our new methodology to calcualte the corresponding UID (is not useful)
    void testingGenomeLevelWithNewMethodology( YRJObject * yrj  ,int differences);

    
    
    map<LONGS, int>  getUnieativeHitsGenus(YRJObject * yrj , int differences );
    
    
    LONGS unieativeLCAKraken(map<LONGS, int> unieativeHits);

    string getLeastCommonLevel(YRJObject * yrj, short krakenShort);
    
    void testingSpeciesLevelWithNewMethodology(YRJObject * yrj   , int differences);
    
    LONGS getSpeciesLevelUID(short shortName);

    LONGS getGenusLevelUID(short shortName);

    string getLeastCommonLevel2(YRJObject * yrj, LONGS finalResultUID);
    
    void calculate_accurcy(string file , unordered_map<string , long> finalResult);
    
    long getElementInTheResult(string element , unordered_map<string , long> finalResult);
    
    
    void testingSpeciesLevelWeightedNewMethodology(YRJObject * yrj   , int differences);
    
    map<LONGS, int>  getUnieativeHitsWeightedSpecies(YRJObject * yrj , int differences );

    
    map<short, int> hits_kmer_with_differences_weighted(LONG kmer , int differences);
    
    void testingSpeciesLevelWithWeightedMethodology(YRJObject * yrj   , int differences);
    void calculate_accurcy_matlab(string file , unordered_map<string , long> finalResult);
    
    void testKrakenOutput(); // this function to test kraken output accuracy
    
    
    
    void testKmerLevelLevel(YRJObject * yrj); // test rach kmer individually
    
    map<  short, short > getHitsandDifferencesKmer(LONG kmer);
    
    
    string getLevel(long lcaUID , long inputUID); // get the level of LCA between two tax IDs
    
    void writeTheFinalResultMap(string s);

    void writeTestKmerLevels();
    
    
    pair<string, string> getMaxMinLevels(YRJObject * yrj ,vector<short> hits);// return the max level and the min level
    
    
    void testKmerLevelLevelMaxMin(YRJObject * yrj);

    

};


#endif /* testingYRJ_hpp */
