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
#include "CoreTaxonomer.hpp"
#include "Configurations.h"
#include "Tree.hpp"

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
    
    
    Tree * tree = new Tree(path_to_the_tree);
    cout << tree->getNumberOfLeaves(0) << endl;

   
    vector<YRJUnieative *>  yrjUnieativeVector  = tree->getYRJUnieariveVector(path_to_the_yrj_databases, pattern);
    
    CoreTaxonomer * core = new CoreTaxonomer(yrjUnieativeVector , pattern);
    
    pair<INT, INT> p =  core->getTheHashedKmer(1);
    
    cout << p.first << endl;
    
    
    
    return 0;
}



