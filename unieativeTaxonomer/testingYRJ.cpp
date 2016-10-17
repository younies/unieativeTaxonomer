//
//  testingYRJ.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/15/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "testingYRJ.hpp"



void Tester::testYRJvector(){
    ifstream simBA5Stream("/export1/project/hondius/testingUnieative/dna.txt");
    ifstream headerSimBA5Stream(this->path_to_simBA5_headers);
    
    
    
    string header;
    string DNA;
    
    
    getline(headerSimBA5Stream, header  );
    getline(simBA5Stream, DNA );
    getline(simBA5Stream, DNA );
    
    
    
    YRJObject yrj(DNA);
    
    
    
    for(auto kmer : yrj.kmersVector)
        cout << kmer << endl;
    
}
