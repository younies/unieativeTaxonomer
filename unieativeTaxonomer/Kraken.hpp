//
//  Kraken.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/3/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef Kraken_hpp
#define Kraken_hpp

#include "headers.h"
#include "YRJObject.hpp"
#include "headers.h"
#include "TreeNode.hpp"
#include "Tree.hpp"
#include "Hash.hpp"
#include "countingTheHistogramDifferences.hpp"


struct shortRead{
    YRJObject * yrj;
    INT UID;
    SHORT shortUID;
};

class Unieative {
    Tree * tree; // the tree that are going to perform the analysis
    vector<Hash *> Hashes;
    vector<pair<YRJObject * , INT> > shortReads; // each short read with the corresponding UID
    
public:
    //Unieative( vector<shortRead> shortReads);
    //Unieative(vector<string>& shortReads );
    
    SHORT getLCA(  LONG kmer , int & differences);
    
    
    vector<INT> getFinalUIDs(YRJObject * yrjObject, int differences);
    
};


#endif /* Kraken_hpp */
