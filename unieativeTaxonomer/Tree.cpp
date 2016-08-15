//
//  Tree.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/11/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "Tree.hpp"

Tree::Tree(string path)
{
    
    ifstream inputStream(path);
    
    
    this->treeNodesVector.resize(this->treeSize);
    string line;
    
    while (getline(inputStream, line))
    {
        stringstream liness(line);
        short shortName;
        bool tagged;
        LONG uid, parentUid;
        
        liness >> uid;
        liness >> parentUid;
        liness >> shortName;
        liness >> tagged;
        
        this->treeNodesVector[shortName]  = new TreeNode( uid ,  parentUid ,  shortName ,  tagged);
        
    }
    
    inputStream.close();
    
    /**
    //set the mapper between the short names and uids
    this->fromShortNameToUid.resize(this->treeNodesVector.size() + 1000 , -1);
    for(LONGS i = 0 , n = this->treeNodesVector.size() ; i < n ; ++i )
        this->fromShortNameToUid[this->treeNodesVector[i]->shortName] = this->treeNodesVector[i]->uid;
    */
}






short Tree::getParentShortName(short shortName)
{
    return this->treeNodesVector[shortName]->parentShortName;
}





