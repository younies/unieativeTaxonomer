//
//  comments.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/21/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "headers.h"

// insert code here...

/*
 Tree * tree = new Tree(path_to_the_tree);
 cout << tree->getNumberOfLeaves(0) << endl;
 
 
 
 vector<YRJObject*> yrjObj = tree->getYRJobjects(path_to_the_yrj_databases);
 
 
 vector<YRJObject*>::const_iterator first = yrjObj.begin() + 0;
 vector<YRJObject*>::const_iterator last = yrjObj.begin() + 10;
 vector<YRJObject*> newVec(first, last);
 
 
 
 CoreTaxonomer * core = new CoreTaxonomer(yrjObj , pattern);
 
 core->writeTheCoreData(path_to_unieative_all_hashed);
 
 HashedNode p =  core->getTheHashedKmer(1);
 
 cout <<sizeof(p) << endl;
 
 
 
 YRJObject * getRandom = new YRJObject(pathe , 0);
 
 getRandom->fillTheKmersVector();
 
 
 ifstream ifsMil(million);
 ifstream ifsTenMil(tenMillion);
 
 vector<LONG> tenMillionsInd;
 vector<LONG> millionsInd;
 
 string temp;
 while (getline(ifsMil, temp))
 {
 millionsInd.push_back( getLong( temp) );
 }
 
 while (getline(ifsTenMil, temp))
 {
 tenMillionsInd.push_back( getLong( temp) );
 }
 
 vector<LONG> tenMillions;
 vector<LONG> millions;
 
 
 for (LONGS i = 0  , n = tenMillionsInd.size() ;  i < n ; ++i)
 {
 tenMillions.push_back( getRandom->kmersVector[  tenMillionsInd[i] ] );
 }
 
 for (LONGS i = 0 , n = millionsInd.size() ; i < n ; ++i)
 {
 millions.push_back( getRandom->kmersVector[  millionsInd[i] ] );
 }
 
 
 cout << millions.size() << endl;
 cout << tenMillions.size() << endl;
 
 sort(tenMillions.begin(), tenMillions.end());
 sort(millions.begin(), millions.end());
 
 
 
 writeFile( millions , 31 , writeIN + "randomMillion.yrj");
 writeFile( tenMillions , 31 , writeIN + "randomTenMillion.yrj");
 
 */



void writeFile(vector<LONG> & kmers , LONG kmerSize , string path)
{
    ofstream output_file(path);
    
    
    LONGS kmerNumbers = kmers.size();
    
    output_file.write( (char *) &kmerSize ,  sizeof(LONG));
    output_file.write( (char *) &kmerNumbers ,  sizeof(LONG));
    
    for (LONGS i = 0 ; i < kmerNumbers ; ++i)
    {
        output_file.write( (char *) &kmers[i] ,  sizeof(LONG));
    }
    
    output_file.flush();
    output_file.close();
    
}


LONGS getLong( string s)
{
    stringstream ss(s);
    LONGS ret;
    ss >> ret;
    return ret;
}


