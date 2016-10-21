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
    
    
    this->treeNodesVector.resize(11000 );
    //this->fromShortNameToUid.resize(this->treeSize , -1);
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
        liness >> parentShortName;
        liness >> tagged;
        
        this->treeNodesVector[shortName]  = new TreeNode( uid ,  parentUid ,  shortName , parentShortName,  tagged);
        
    }
    
    inputStream.close();
    
    //connect parents
    this->connectChildren();
    
    //build the leaves numbers
    //this->buildTheNumberOfLeaves();
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
    cout << "short n  " << shortName << endl;
    cout << this->treeNodesVector.size() << endl;
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





short Tree::getLCA_for_child_node( TreeNode * node , vector<short> & hitted_nodes )
{
    if( binary_search(hitted_nodes.begin(), hitted_nodes.end(), node->shortName))
        return node->shortName;
    if( this->treeNodesVector[node->shortName]->children.size() == 0)
        return -1;
    vector<short> LCAs;
    
    for(TreeNode * child : node->children)
    {
        short lca = getLCA_for_child_node(child , hitted_nodes ) ;
        if( lca != -1)
        {
            LCAs.push_back(lca);
        }
    }
    
    switch (LCAs.size()) {
        case 0:
            return -1;
            break;
        case 1:
            return LCAs[0];
            break;
            
        default:
            return node->shortName;
            break;
    }
    
    return  -1;
}





G_Statistics Tree::calculateG_Statistics(LONG kmer , short nodeShortName , vector<short> & hitted_nodes  )
{
    G_Statistics gStat;
    
    gStat.kmer = kmer;
    
    gStat.number_of_hitted_leaves = this->numberOfLeaves[nodeShortName];
    
    gStat.demoneratorGX = 0;
    
    for( TreeNode * child : this->treeNodesVector[nodeShortName]->children)
    {
        short lca = getLCA_for_child_node(child , hitted_nodes);
        
        if(lca >= 0)
            gStat.demoneratorGX += this->numberOfLeaves[lca];
    }
    
    gStat.number_of_hitted_leaves = hitted_nodes.size();
    
    gStat.LCA_global_shortName = nodeShortName;
    
    gStat.LCA_global_Uid = this->treeNodesVector[nodeShortName]->uid;
    
    gStat.GX = (double)gStat.number_of_leaves / (double)gStat.demoneratorGX;
 
    return gStat;
}




short Tree::getTowLCA(short first , short second)
{
    while(first != second)
    {
        if(first > second)
            first = this->getParentShortName(first);
        else
            second = this->getParentShortName(second);
    }
    
    return first;
}


short Tree::getGlobalLCA( vector<short> &  hitted_nodes)
{
    LONGS vecSize = hitted_nodes.size();
    
    if(vecSize == 0)
    {
        cout << "error in the get Global LCA 22\n";
        return  -1;
    }
    else if(vecSize == 1)
        return hitted_nodes[0]; // return the index of that element\
    
    short ret =  hitted_nodes[0];
    
    for(LONGS i = 1 ; i < vecSize  ; ++i )
        ret = getTowLCA(ret, hitted_nodes[i]);
    
    return ret;

    
}



short Tree::getGlobalLCA( vector<pair< short , short> > & vectorOfResults)
{
    LONGS vecSize = vectorOfResults.size();
    
    if(vecSize == 0)
    {
        cout << "error in the get Global LCA\n";
        return  -1;
    }
    else if(vecSize == 1)
        return vectorOfResults[0].first; // return the index of that element
    
    
    short ret =  vectorOfResults[0].first;
    for(LONGS i = 1 ; i < vecSize  ; ++i )
        ret = getTowLCA(ret, vectorOfResults[i].first);
    
    return ret;
    
}



LONGS Tree::getNumberOfLeaves(short shortName)
{
    return this->numberOfLeaves[shortName];
}


vector<YRJObject *> Tree::getYRJobjects(string path )
{
    vector<YRJObject *> ret;
    
    for ( TreeNode * node : this->treeNodesVector)
    {
        if(node != NULL && node->tagged)
        {
            ostringstream ss;
            ss << node->uid;
            ret.push_back( (new YRJObject(path + ss.str() + ".yrj" , node->shortName)) );
        }
    }
    
    return ret;
}







LONGS Tree::getTheUIDFromShort(short shortName)
{
    return this->treeNodesVector[shortName]->uid;
}







