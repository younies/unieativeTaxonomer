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

class Tester {
    vector<YRJObject* > yrjVector;
    
    const string path_to_simBA5 = "/export1/project/hondius/testingUnieative/accuracy/simBA5_accuracy.fa";
    const string path_to_simBA5_headers = "/export1/project/hondius/testingUnieative/accuracy/mappingGI.txt";
    
public:
    
    
    void testYRJvector();
};


#endif /* testingYRJ_hpp */
