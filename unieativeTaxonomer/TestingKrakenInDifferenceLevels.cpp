//
//  TestingKrakenInDifferenceLevels.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/9/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "TestingKrakenInDifferenceLevels.hpp"










//implementing the sigmBA5 testing file
void testEntiresimBA5(string path_to_simBA5  , string path_to_simBA5_headers, string path_to_the_names_dmp_file , string path_to_the_nodes_dmp_file , string path_Gi, string path_to_the_output_file , string path_to_hash_database , vector<string> hashesPattern , string pruined_tree , int numOfDifferences)
{
    ofstream * resultFile = new ofstream(path_to_the_output_file); // creating the output file reference
    
    vector<YRJObject *> yrjObjects; // the yrj objects references in a vector
    
    ifstream * simBA5_file = new ifstream(path_to_simBA5);
    ifstream * simBA5_headers_file = new ifstream(path_to_simBA5_headers);
    
    string fastaDNA;
    string fastaHeader;
    
    
    
    BigTree * bigTree = new BigTree(path_to_the_names_dmp_file ,path_to_the_nodes_dmp_file, path_Gi );
    
    Tree * tree = new Tree(pruined_tree);
    
    
    Unieative * unieative = new Unieative(tree , hashesPattern , path_to_hash_database);
    
    
    
    while (getline(*simBA5_headers_file, fastaHeader)) {
        getline(*simBA5_file, fastaDNA);
        getline(*simBA5_file, fastaDNA);
        
        
        YRJObject * yrjObject = new YRJObject(fastaDNA);
        
        yrjObject->uid = bigTree->getUIDFromFastaHeaderGI(fastaHeader);
        yrjObjects.emplace_back(yrjObject);
        
    }
    
    vector<LONGS> globalResults(numOfDifferences+1 , 0);
    
    
    for(YRJObject*  yrj : yrjObjects){

        testOutYRJfile(resultFile, tree, bigTree, yrj, unieative, numOfDifferences, globalResults);
    }
    
    
}








void testOutYRJfile( ofstream * writeFile , Tree * tree, BigTree * bigTree, YRJObject * yrjObject  ,  Unieative * unieative ,int numOfDifferences , vector<LONGS> &globalResult )
{
    
    
    for (int  currDiff = 0 ; currDiff <= numOfDifferences ; ++currDiff) {
        int testCorrecteness = 0;
        LONGS assignedUid = unieative->getFinalUIDs(yrjObject, currDiff);
        
        if(bigTree->isBothInGenusLevel(assignedUid, yrjObject->uid))
            testCorrecteness = 1;
        
        *writeFile << testCorrecteness << "\t";
        globalResult[currDiff] += testCorrecteness;
    }
    
    *writeFile << endl;
    
    for(auto num : globalResult)
        *writeFile << num << "\t";
    *writeFile << endl;
    
}



