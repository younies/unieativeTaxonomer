//
//  Count_differences_support.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 9/19/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "countingTheHistogramDifferences.hpp"


vector<vector< LONGS> >  addLocalToGlobal( vector<vector<LONGS> > &globalCounter , const vector<vector<vector<SHORT>>> & localCounter)
{
    vector<vector< LONGS> > ret = globalCounter;
    for (int i = 0 , n = (int)globalCounter.size() ; i < n ; ++i)
    {
        for (int j = 0 ,  m = (int) globalCounter[i].size() ; j < m ; ++j)
        {
            ret[i][j] = localCounter[i][j].size();
            globalCounter[i][j] += ret[i][j];
        }
    }
    
    
    return ret;
    
    
}



void writeMatrixInFile(ofstream * writtingFile , vector<vector<LONGS> > data)
{
    for (int i = 0 ,   n = (int)data.size() ;  i <  n ;  ++ i )
    {
        for (int j = 0 ,   m = (int)data[i].size();   j <  m  ;  ++j ) {
            
            (*writtingFile) << data[i][j] << "\t";
        }
        *writtingFile << endl;
    }
}



