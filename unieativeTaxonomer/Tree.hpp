//
//  Tree.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 8/11/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef Tree_hpp
#define Tree_hpp
#include "YRJObject.hpp"
#include "headers.h"
#include "TreeNode.hpp"
#include "BigTree.hpp"

/**
 
 this class is the abstraction of the pruined tree (the tree after applying pruining on it), and 
 each node in the tree is represented by class TreeNode.
 
 Most of the operations on this class are optimized because of the hierarchical of the pruining tree.
 */

class Tree
{
    vector< TreeNode * > treeNodesVector;
    //vector<LONGS> fromShortNameToUid;
    const long treeSize = 11000l;
    vector<LONGS> numberOfLeaves;
    void connectChildren();
    LONGS setNumberOfLeaves(TreeNode * node);
    const string species    = "species";
    const string genus      = "genus";
    
    //vector for the genus level and species  level
    
    

public:
    BigTree * bigTree;

    short root = 0; /// should always hold the short for the root
    
    /*
     Constructor and it needs the path to the pruined tree.
     */
    Tree(string pathToTheTree);
    
    //return the shortUID for the parent by taking the shortUID for a node.
    short getParentShortName(short shortName);
    
    //This method fill the number of leaves vector which each index is the shortUID
    //and the value is the number of leaves.
    void buildTheNumberOfLeaves( );
    
    //to get the g_statistics domenrator
    
    // return the lowest common ancesstor between two short UIDs
    short getTowLCA(short first , short second);
    
    // return the lowest common ancesstor between vector of short UIDs
    short getGlobalLCA( vector<short> &  hitted_nodes);
    
    // take vector of pairs of indices and differences and return Global LCA
    short getGlobalLCA( vector<pair< short , short> > & vectorOfResults);

    //return the number of leaves under this short node
    LONGS getNumberOfLeaves(short shortName);
    
    //return the number of leaves
    short getLCA_for_child_node( TreeNode * node , vector<short> & hitted_nodes );
    
    
    // This method return the G statistics analysis in an object
    G_Statistics calculateG_Statistics(LONG kmer , short nodeShortName  , vector<short> & hitted_nodes );
    
    
    //retun the actual UID from the short UID
    LONGS getTheUIDFromShort(short shortName);
    
    
    //unused!!
    vector<YRJObject *> getYRJobjects(string path );
    
    // this function return the speciest parent short number
    short getSpeciesParent(YRJObject * yrj);

    // this function return the Genus parent short number
    short getGenusParent(YRJObject * yrj);
    
    // this function return the speciest parent short number
    short getSpeciesParent(short shortName);
    
    // this function return the Genus parent short number
    short getGenusParent(short shortName);
    
    // take the binary tree with the short UIDs and the number of hits and return the Kraken final result.
    short getTheMaximumKRAKENhit(map<short, int> & hitNumbers);
    
    // unused!!
    void connectSpeciesGenusLevels(BigTree * bigTree);
    
};
#endif /* Tree_hpp */
