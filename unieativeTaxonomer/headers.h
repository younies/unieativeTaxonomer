//
//  headers.h
//  building_vectors
//
//  Created by Younies Mahmoud on 7/15/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef headers_h
#define headers_h

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <err.h>
#include <errno.h>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdint.h>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sysexits.h>
#include <unistd.h>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <stack>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <memory.h>
#include <fstream>

//headaers for the program itself




using namespace std;

typedef unsigned long LONG;
typedef long LONGS;




//abstracting the new hashed kmer
struct HashedNode{
    short index;
    unsigned int rawKmer; //the non hashed part in the kmer
    unsigned int hashedKmer; // the part that is hashed in the kmer
    
    bool operator< ( const HashedNode& y) {
        return std::tie(this->rawKmer, this->index , this->hashedKmer ) < std::tie(y.rawKmer, y.index , y.hashedKmer);
    }
};

bool hashedNodeCompare(HashedNode &lhs, HashedNode &rhs) { return lhs < rhs; } // for the binartysearch or the sorting algorithms
//end of the hasehd kmer implementation


struct Node
{
    LONG uid;
    LONG parentUID;
    LONGS parentIndex;
    LONGS myselfIndex;
    vector<LONG> children;
};

struct G_Statistics
{
    LONG kmer;
    LONGS number_of_leaves;
    LONGS number_of_hitted_leaves;
    Node LCA_global;
    LONGS demoneratorGX;
    double GX;
    
};

struct Uid_Value {
    LONGS uid;
    LONGS value;
};

#endif /* headers_h */
