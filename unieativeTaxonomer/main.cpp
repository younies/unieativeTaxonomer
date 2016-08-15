//
//  main.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "headers.h"
#include "YRJObject.hpp"
#include "YRJUnieative.hpp"
#include "helpers.hpp"

string pattern = "##-#--###---#-#-#-#-#--#-##--##";
string pathe  = "/Users/youniesmahmoud/Desktop/903893.yrj";// = "/export1/project/hondius/newKrakenResearch/databases/kmerDatabase_new_31/all.yrj";

int main(int argc, const char * argv[])
{
    // insert code here...
    /**
    YRJUnieative *node =new YRJUnieative(pathe , pattern);
    
    node->fillTheHashedNodesVector();
    delete node;
    
     */
    
    cout << reverseKmer(1);
    
    return 0;
}



