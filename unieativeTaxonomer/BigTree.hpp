//
//  BigTree.hpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/8/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#ifndef BigTree_hpp
#define BigTree_hpp



#include "headers.h"


//the class that contains all the required data structures for the tree



class BigTree {
    
    //pathes to the node and names files
    string path_nodes;
    string path_names;
    
    //for checking for specific names
    bool is_specific_leaves = false;
    string path_to_specific_names;
    vector<LONGS> specific_nodes_indices_sorted;
    
    vector< LONG> ids;//holds the uid for each node index
    vector<vector<string > > names;//it associate all the names for a node in the tree
    vector<string  > levels; // this specify the level for each node (if it is species, genus .... and so on)
    vector< Node > trie;// holds all the trie nodes starting with the root
    map< LONGS , LONGS > GiMap;
    
    vector<Node> sorted_leafs_df; // this vector holds the indices  of sorted leafs in the depth first order
    LONGS strTolong(string s); // this function convert the string to integer type (long long)
    
    // for constructing the tree
    void construct_ids_names(); // this method construct the ids and names vector
    void constructing_trie(); // this method construct the trie
    void construct_sorted_leafs(LONG index  = 0 ); //this method construct the sorted leafs vector
    
    
public:
    BigTree(string path_names , string path_nodes);// the default constructor
    BigTree(string path_names , string path_nodes , string path_Gi);
    
    
    
    string getLevel(LONGS position); // method for getting the type of an index node
    
    LONGS uid_to_index(LONG uid);// methid that convert from node id to its index
    
    
    vector<LONGS> get_All_Parents(Node node); // return parents indices
    LONGS get_LCA_between_Two_Nodes(Node node1 , Node node2);
    Node get_Global_LCA(vector<Node>  nodes);
    
    
    LONGS get_Number_Of_Node_Leaves(Node node);
    
    LONGS getNumberOfHits(Node node , vector<LONGS> & sorted_indicies);// the indicies should be sorted
    
    void get_hitted_nodes(Node node , vector<LONGS> & sorted_indicies , vector<Node> & hitted_Nodes ); // return list of node hitted under the given node
    
    Node getNodeFromIndex(LONGS index);
    
    
    
    vector<LONG> getSortedLeafsUIDs();
    
    void build_specific_nodes(string path_to_nodes);
    void update_specific_nodes(vector<LONG> UIDs);
    
    bool is_this_Index_node_specified(LONGS index);
    
    
    LONGS get_demonrator_GX(Node node ,vector<LONGS> & indices);
    
    vector<LONGS> get_global_kraken_uids(vector<Node> nodes);
    
    void calculate_index_all_uid_values( vector<Uid_Value> & all_index_values );
    
    LONGS get_vlue_for_parent(vector<Uid_Value> & all_index_values , LONGS index);
    
    string get_level(Node node);
    
    LONGS getUIDFromFastaHeaderGI(string fastaHeader);// return the UID for header such as "simBA.000000308 160901491"
    
    bool isBothInGenusLevel(LONGS uid_first , LONGS uid_second);
    
    LONGS getGenusUID(LONGS uid);
    
    LONGS getTheSpeciesUID(LONGS uid);// this method return the species UID for the taxonomer or -1 if it is upove this UID
    
    string getNextNotNoLevellevel(Node node);
    
};







#endif /* BigTree_hpp */
