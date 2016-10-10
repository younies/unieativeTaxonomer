//
//  BigTreeNewFunctionality.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/9/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "BigTree.hpp"



LONGS BigTree::getUIDFromFastaHeaderGI(string fastaHeader)
{
    stringstream fastaStream(fastaHeader);
    
    LONGS ret;
    
    getline(fastaStream, fastaHeader, ' ');
    
    getline(fastaStream, fastaHeader);
    
    ret = stol(fastaHeader);
    
    return ret;
}





BigTree::BigTree(string path_names , string path_nodes , string path_Gi)
{
    BigTree( path_names ,  path_nodes );
    
    ifstream Gi_FileStream(path_Gi);
    
    string giString;
    string uidString;
    while (getline(Gi_FileStream, giString, ' '))
    {
        getline(Gi_FileStream, uidString);
        
        LONGS gi = stol(giString);
        LONGS uid = stol(uidString);
        
        this->GiMap[gi] = uid;
    }

}
