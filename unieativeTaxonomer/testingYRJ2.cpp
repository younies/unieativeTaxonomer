//
//  testingYRJ2.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/21/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "testingYRJ.hpp"



void Tester::testingSpeciesLevelWithNewMethodology(YRJObject * yrj   , int differences)
{
    auto  unieativeHits = this->getUnieativeHitsSpecies(yrj ,  differences );
    if(unieativeHits.empty()){
        this->finalResult[this->notConsidered]++;
        return;
    }

    auto finalUnieative = this->unieativeLCAKraken(unieativeHits);

    
    if(bigTree->getGenusUID(finalUnieative) < 0)
    {
        finalResult[notNeeded]++;
        return;
    }
    
    auto level = this->getLeastCommonLevel2(yrj, finalUnieative);

    if(this->finalResult.count(level))
        finalResult[level] ++;
    else
        finalResult[level] = 1;
    
    

}



void Tester::testingGenomeLevelWithNewMethodology(YRJObject * yrj   , int differences)
{
    //using my new methodology to make it better
    
    auto  unieativeHits = this->getUnieativeHitsGenus(yrj ,  differences );
    
    if(unieativeHits.empty()){
        this->finalResult[this->notConsidered]++;
        return;
    }
    
    auto finalUnieative = this->unieativeLCAKraken(unieativeHits);
    
    auto level = this->getLeastCommonLevel(yrj, finalUnieative);
    
    if(bigTree->getGenusUID(finalUnieative) < 0)
    {
        finalResult[notNeeded]++;
        return;
    }

    
    
    if(bigTree->getGenusUID(finalUnieative) < 0)
    {
        finalResult[notNeeded]++;
        return;
    }
    
    
    if(this->finalResult.count(level))
        finalResult[level] ++;
    else{
        finalResult[level] = 1;
    }
    //cout << level << endl;
    /*
     
    string s;
    s.push_back( (char)(differences + '0'));

    ofstream * result = new ofstream(this->result + s + ".out" );
    for(auto res : this->finalResult)
    {
        *result << res.first << "\t" << res.second << endl;
    }
    */
}


//Kraken
void Tester::testingGenomeLevel(YRJObject * yrj   , int differences)
{
    /*
    if(!this->isKrakenCatch(yrj))
    {
        this->finalResult[this->notConsidered] ++;
        return;
    }
     */
    auto hitNumbers = this->getKrakenLCAs(yrj, differences);
    
    //to get kraken Final Taxonomy
    
    if(hitNumbers.size() == 0)
    {
        this->finalResult[this->notConsidered] ++;
        return;
    }
    
    auto krakenShort = this->pruinedTree->getTheMaximumKRAKENhit(hitNumbers);
    
    
    auto UID = this->pruinedTree->getTheUIDFromShort(krakenShort);
    
    auto genusUID = this->bigTree->getGenusUID(UID);
    
    if(genusUID == -1)
    {
        finalResult[this->notNeeded]++;
        return;
    }
    
    auto level = this->getLeastCommonLevel(yrj, krakenShort);
    
    if(this->finalResult.count(level))
        finalResult[level] ++;
    else
        finalResult[level] = 1;
    
}




bool Tester::isKrakenCatch(YRJObject * yrj )
{
    //to find all the LCAs
    for(auto kmer: yrj->kmersVector)
    {
        auto hits = this->hits_kmer_with_differences(kmer, 0);
        
        if(hits.size() != 0)
            return true;
    }
    
    return false;

}

map<short, int>  Tester::getKrakenLCAs(YRJObject * yrj , int differences)
{
    map<short, int> hitNumbers;
    //to find all the LCAs
    for(auto kmer: yrj->kmersVector)
    {
        auto hits = this->hits_kmer_with_differences(kmer, differences);
        
        if(hits.size() == 0){
            continue;
        }
        
        if(hits.size() > 0)
        {
            auto lca = this->pruinedTree->getGlobalLCA(hits);
            
            
            
            
            
            
            if(hitNumbers.count(lca))
                hitNumbers[lca]++;
            else
                hitNumbers[lca] = 1;
        }
    }
    return hitNumbers;
}


map<LONGS, int>  Tester::getUnieativeHitsSpecies(YRJObject * yrj , int differences )
{
    map<LONGS, int> unieativeHits;
    
    for(auto kmer : yrj->kmersVector)
    {
        vector<short> hits = this->hits_kmer_with_differences(kmer, differences);
        
        if(hits.size() == 0){
            continue;
        }
        
        set<LONGS> speciesHits;
        
        for (auto hit: hits)
        {
           // cout << "start \n";
            //auto hit = this->pruinedTree->getTheUIDFromShort(hito);
            
          //  cout << "end \n";
            if(this->getSpeciesLevelUID( hit) != -1)
                speciesHits.insert( this->getSpeciesLevelUID( hit));
            else
                cerr << "big problem 44 \n" << hit <<endl;
        }
        
        for(auto species: speciesHits)
        {
            if(unieativeHits.count(species))
                unieativeHits[species]++;
            else
                unieativeHits[species] = 1;
        }
        
    }
    
    
    return unieativeHits;

}



map<LONGS, int>  Tester::getUnieativeHitsGenus(YRJObject * yrj , int differences )
{
    map<LONGS, int> unieativeHits;
    
    
    
    for(auto kmer : yrj->kmersVector)
    {
        auto hits = this->hits_kmer_with_differences(kmer, differences);
        
        
        if(hits.size() == 0){
            continue;
        }
        
        set<LONGS> genusHits;
        
        for (auto hit: hits)
        {
            if(hit < 0)
            {
                cout << "weired error  " << hit << "\n";
                break;
            }
            auto genusUID = getGenusLevelUID(hit);
            if(genusUID != -1)
                genusHits.insert( genusUID);
            else
                cerr << "big problem \n" << hit <<endl;
        }
        
        
        for(auto genusHit : genusHits)
        {
            if(unieativeHits.count(genusHit))
                unieativeHits[genusHit]++;
            else
                unieativeHits[genusHit] = 1;
        }
        
    }
    
    return unieativeHits;
     
}



LONGS Tester::unieativeLCAKraken(map<LONGS, int> unieativeHits)
{
    int maximum = 0;
    
    vector<Node> maxNodes;
    for(auto hit: unieativeHits)
        maximum = max(maximum , hit.second);
    
    for(auto unieative: unieativeHits)
        if(unieative.second == maximum)
            maxNodes.emplace_back( this->bigTree->getNodeFromIndex( this->bigTree->uid_to_index(  unieative.first)));
    
    
    auto finalNode = this->bigTree->get_Global_LCA(maxNodes);
    
    return finalNode.uid;

    
}


string Tester::getLeastCommonLevel2(YRJObject * yrj, LONGS finalResultUID)
{
    auto krakenNode  = this->bigTree->getNodeFromIndex(this->bigTree->uid_to_index(finalResultUID));
    auto indexYRJ = this->bigTree->uid_to_index(yrj->uid);
    
    auto yrjNode     = this->bigTree->getNodeFromIndex(indexYRJ);
    
    auto testingLevelIndex = this->bigTree->get_LCA_between_Two_Nodes(krakenNode, yrjNode);
    
    auto levelNode = this->bigTree->getNodeFromIndex(testingLevelIndex);
    
    string level = this->bigTree->get_level(levelNode);
    
    if(level == "no")
    {
        level = this->bigTree->getNextNotNoLevellevel(levelNode);
    }
    
    //if( bigTree->getGenusUID(levelNode.uid)  < 0 )
      //  return this->notNeeded;
    
    return level;
}



string Tester::getLeastCommonLevel(YRJObject * yrj, short krakenShort)
{
    auto krakenUID   = this->pruinedTree->getTheUIDFromShort(krakenShort);
    
    auto krakenNode  = this->bigTree->getNodeFromIndex(this->bigTree->uid_to_index(krakenUID));
    
    auto indexYRJ = this->bigTree->uid_to_index(yrj->uid);
    
    auto yrjNode     = this->bigTree->getNodeFromIndex(indexYRJ);
    
    auto testingLevelIndex = this->bigTree->get_LCA_between_Two_Nodes(krakenNode, yrjNode);
    
    auto levelUID = this->bigTree->getNodeFromIndex(testingLevelIndex);
    
    string level = this->bigTree->get_level(levelUID);
    
    if(level == "no")
        level = this->bigTree->getNextNotNoLevellevel(levelUID);
    
    return level;

}






LONGS Tester::getSpeciesLevelUID(short shortName)
{
    auto uid = this->pruinedTree->getTheUIDFromShort(shortName);
    
    return this->bigTree->getTheSpeciesUID(uid);
}


LONGS Tester::getGenusLevelUID(short shortName)
{
    auto uid = this->pruinedTree->getTheUIDFromShort(shortName);

    return this->bigTree->getGenusUID(uid);
}



void Tester::calculate_accurcy(string file)
{
    ofstream output(file);
    
    
    double total = 0.0;
    
    for(auto element : finalResult)
    {
        
        if(element.first == notConsidered)
            continue;
        if(element.first == notNeeded)
            continue;
        
        total+= element.second;
    }
    
    
    double speciesAccuracy = 0.0;
    
    
    speciesAccuracy += getElementInTheResult("species");
    speciesAccuracy += getElementInTheResult("subspecies");
    
    
    
    
    double genusAccuracy = speciesAccuracy ;
    genusAccuracy +=getElementInTheResult ("genus");
    
    
    double not_counted = getElementInTheResult(notConsidered) + getElementInTheResult(notNeeded);
    
    
    output << "species:\t" << speciesAccuracy/total << endl;
    output << "genus:\t" << genusAccuracy/total << endl;
    output << "not used\t" << not_counted/(total + not_counted) << endl;
    
    output.close();
    
    
    
    
}



long Tester::getElementInTheResult(string element)
{
    if(finalResult.count(element) )
        return finalResult[element] ;

    return 0;
}


