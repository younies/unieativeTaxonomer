//
//  BigTreeExtend.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/10/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//



#include "BigTree.hpp"

bool BigTree::isBothInGenusLevel(LONGS uid_first , LONGS uid_second){

    auto firstGenus = getGenusUID(uid_first);
    if(firstGenus == -1)
        return false;
    
    auto secondGenus = getGenusUID(uid_second);
    if(secondGenus == firstGenus)
        return true;
    else
        return false;

}




LONGS BigTree::getGenusUID(LONGS uid){
    
    if(uid < 2)
        return -1;
    
    while (uid != 1) {
        auto index = this->uid_to_index(uid);
        auto node  = this->getNodeFromIndex(index);
        auto level = get_level(node);
        if(level == "genus")// we should check this
            return uid;
        uid = node.parentUID;
    }
    
    return -1;
    
}








LONGS BigTree::getTheSpeciesUID(LONGS uid)
{
    LONGS index = this->uid_to_index(uid);
    
    
    while (index)
    {
        if(levels[index] == "species")
            return this->trie[index].uid;
        
        index = trie[index].parentIndex;
    }
    
    return -1;
    
}

string BigTree::getLevel(LONGS position)
{
    return levels[position];
}
