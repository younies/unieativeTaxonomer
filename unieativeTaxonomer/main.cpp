//
//  main.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "headers.h"
#include "YRJObject.hpp"
#include "helpers.hpp"
#include "CoreTaxonomer.hpp"
#include "Configurations.h"
#include "Tree.hpp"

string pattern = "##-#--###---#-#-#-#-#--#-##--##";
//string pathe  = "/Users/youniesmahmoud/Desktop/903893.yrj";// = "/export1/project/hondius/newKrakenResearch/databases/kmerDatabase_new_31/all.yrj";

int main(int argc, const char * argv[])
{
    // insert code here...
    
    
   Tree * tree = new Tree(path_to_the_tree);
    cout << tree->getNumberOfLeaves(0) << endl;

   
    
    vector<YRJObject*> yrjObj = tree->getYRJobjects(path_to_the_yrj_databases);
    
    /*
    vector<YRJObject*>::const_iterator first = yrjObj.begin() + 0;
    vector<YRJObject*>::const_iterator last = yrjObj.begin() + 10;
    vector<YRJObject*> newVec(first, last);
    */
    CoreTaxonomer * core = new CoreTaxonomer(yrjObj , pattern);
    
    HashedNode p =  core->getTheHashedKmer(1);
    
    cout <<sizeof(p) << endl;
    
    
    
    return 0;
}



