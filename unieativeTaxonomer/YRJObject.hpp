//
//  YRJObject.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef YRJObject_hpp
#define YRJObject_hpp

#include "headers.h"


class YRJObject
{
    LONG kmerLength;
    LONG numOfKmers;
    string path_to_file;
    short index;
    
public:
    YRJObject(string path);
    ~YRJObject();
    
    vector<LONG>  kmersVector;
    short getIndex();
    LONGS getNumOfKmers();
    
    void fillTheKmersVector();
    
    void  clearTheCompleteKmers();
};

#endif /* YRJObject_hpp */
