//
//  main.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "headers.h"
#include "YRJObject.hpp"


string pathe = "/Users/youniesmahmoud/Desktop/903893.yrj";

int main(int argc, const char * argv[]) {
    // insert code here...
    
    YRJObject *yrj = new YRJObject(pathe);
    
    delete yrj;
    
    
    cerr << "happend\n";
    return 0;
}
