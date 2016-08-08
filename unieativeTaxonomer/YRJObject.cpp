//
//  YRJObject.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "YRJObject.hpp"


YRJObject::YRJObject(string path)
{
    this->path_to_file = path_to_file;
    ifstream fileStream(path_to_file);
    if(!fileStream.good())
    {
        cerr<< "file not found!!!!\n" + path;
        return;
    }
    fileStream.read( (char *)&this->kmerLength  , sizeof(LONG));
    fileStream.read( (char *) &this->numOfKmers , sizeof(LONG));
    this->kmersVector.resize(this->numOfKmers);
    
    fileStream.read((char *)kmersVector[0], this->numOfKmers * sizeof(LONG));
    fileStream.close();
}