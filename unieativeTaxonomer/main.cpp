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

string pattern = "##-#--###---#-#-#-#-#--#-##--##";
string pathe  = "/export1/project/hondius/newKrakenResearch/databases/kmerDatabase_new_31/all.yrj";

int main(int argc, const char * argv[])
{
    // insert code here...
    
    YRJUnieative *node =new YRJUnieative(pathe , pattern);
    
    delete node;
    
    return 0;
}



