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
    LONGS kmerLength;
    LONGS numOfKmers;
    string path_to_file;
    short index;
    
    ifstream * fileStream;
public:
    
    
    
    YRJObject(string path , short index);
    ~YRJObject();
    
    vector<LONG>  kmersVector;
    short getIndex();
    LONGS getNumOfKmers();
    
    void fillTheKmersVector();
    
    void  clearTheCompleteKmers();
    
    bool openFileStream();
    LONG readAKmer();
    void closeFileStream();
    
};

#endif /* YRJObject_hpp */
