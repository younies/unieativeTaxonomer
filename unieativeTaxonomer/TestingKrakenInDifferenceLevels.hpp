//
//  TestingKrakenInDifferenceLevels.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/9/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef TestingKrakenInDifferenceLevels_hpp
#define TestingKrakenInDifferenceLevels_hpp

#include "BigTree.hpp"
#include "headers.h"
#include "YRJObject.hpp"
#include "Tree.hpp"
#include "Hash.hpp"
#include "HashedNode.hpp"
#include "Kraken.hpp"



void testEntiresimBA5(string path_to_simBA5  , string path_to_simBA5_headers, string path_to_the_names_dmp_file , string path_to_the_nodes_dmp_file , string path_Gi, string path_to_the_output_file );


void testOutYRJfile( ofstream * writeFile , vector<Hash*> &hashes, Tree * tree, BigTree * bigTree, YRJObject * yrjObject  ,  Unieative * unieative ,int numOfDifferences );


#endif /* TestingKrakenInDifferenceLevels_hpp */
