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
    
    
    this->treeNodesVector.resize(this->treeSize , NULL);
    string line;
    
    while (getline(inputStream, line))
    {
        stringstream liness(line);
        short shortName;
        short parentShortName;
        bool tagged;
        LONG uid, parentUid;
        
        liness >> uid;
        liness >> parentUid;
        liness >> shortName;
        liness >> tagged;
        
        this->treeNodesVector[shortName]  = new TreeNode( uid ,  parentUid ,  shortName , parentShortName,  tagged);
        
    }
    
    inputStream.close();
    
    //connect parents
    this->connectChildren();
}


void Tree::connectChildren()
{
    for (LONGS i = 0 , n = this->treeNodesVector.size(); i < n ; ++i)
    {
        if(this->treeNodesVector[i] != NULL && this->treeNodesVector[i]->parentShortName != this->treeNodesVector[i]->shortName)//to avoid the loop in the root
            this->treeNodesVector[ this->treeNodesVector[i]->parentShortName ]->children.push_back(this->treeNodesVector[i]);
    }
}





short Tree::getParentShortName(short shortName)
{
    return this->treeNodesVector[shortName]->parentShortName;
}



LONGS Tree::setNumberOfLeaves(TreeNode * node)
{
    LONGS numOfLeaves = 0;
    if(node->tagged)
        numOfLeaves += 1;
    
    for (LONGS i = 0 , n = node->children.size(); i < n ; ++i)
    {
        numOfLeaves += setNumberOfLeaves(node->children[i]);
    }
    
    this->numberOfLeaves[node->shortName] = numOfLeaves;
    
    return numOfLeaves;
    
}



void Tree::buildTheNumberOfLeaves( )
{
    this->numberOfLeaves.resize(this->treeNodesVector.size() , 0);
    
    setNumberOfLeaves(this->treeNodesVector[0]);
}































