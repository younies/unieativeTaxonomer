//
//  YRJObject.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "YRJObject.hpp"


YRJObject::YRJObject(string path ,short index)
{
    this->index = index;
    this->path_to_file = path;
    ifstream fileStream(path_to_file);
    if(!fileStream.is_open())
    {
        cerr<< "file not found!!!!\n" + path;
        return;
    }
    fileStream.read( (char *)&this->kmerLength  , sizeof(LONG));
    fileStream.read( (char *) &this->numOfKmers , sizeof(LONG));
    fileStream.close();
}


YRJObject::~YRJObject()
{
    this->kmersVector.clear();
    cout << "YRJ object is destroyed\n";
}





short YRJObject::getIndex()
{
    return this->index;
}



LONGS YRJObject::getNumOfKmers()
{
    return this->numOfKmers;
}


void YRJObject::clearTheCompleteKmers()
{
    this->kmersVector.clear();
}



void YRJObject::fillTheKmersVector()
{
    clearTheCompleteKmers();
    ifstream fileStream(path_to_file);
    if(!fileStream.is_open())
    {
        cerr<< "file not found!!!! from filling \n" + this->path_to_file;
        return;
    }
    fileStream.read( (char *)&this->kmerLength  , sizeof(LONG));
    fileStream.read( (char *) &this->numOfKmers , sizeof(LONG));
    this->kmersVector.resize(this->numOfKmers);
    
    fileStream.read((char *)&kmersVector[0], this->numOfKmers * sizeof(LONG));
    fileStream.close();

}




bool YRJObject::openFileStream()
{
    ifstream fileStream(path_to_file);
    this->fileStream = &fileStream;
    if(!this->fileStream->is_open())
    {
        cerr<< "file not found!!!! from filling \n" + this->path_to_file;
        return false;
    }
    this->fileStream->read( (char *)&this->kmerLength  , sizeof(LONG));
    this->fileStream->read( (char *) &this->numOfKmers , sizeof(LONG));
    
    return true;
}




LONG YRJObject::readAKmer()
{
    LONG ret;
    this->fileStream->read( (char *)&ret  , sizeof(LONG));
    return ret;
}


void YRJObject::closeFileStream()
{
    this->fileStream->close();
}


