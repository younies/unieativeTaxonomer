//
//  testingYRJ2.cpp
//  unieativeTaxonomer
//
//  Created by Younies Mahmoud on 10/21/16.
//  Copyright © 2016 Younies Mahmoud. All rights reserved.
//

#include "testingYRJ.hpp"

void Tester::testingSpeciesLevelWithWeightedMethodology(YRJObject * yrj   , int differences)
{
    auto  unieativeHits = this->getUnieativeHitsWeightedSpecies(yrj ,  differences );

    if(unieativeHits.empty()){
        this->finalResult[this->notConsidered]++;
        return;
    }
    
    auto finalUnieative = this->unieativeLCAKraken(unieativeHits);

    cout << finalUnieative << endl;
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

void Tester::testingSpeciesLevelWithNewMethodology(YRJObject * yrj   , int differences)
{
    auto  unieativeHits = this->getUnieativeHitsSpecies(yrj ,  differences );
    if(unieativeHits.empty()){
        this->finalResult[this->notConsidered]++;
        return;
    }

    auto finalUnieative = this->unieativeLCAKraken(unieativeHits);

    cout << finalUnieative << endl;
    
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
    long sum = 0;
    for(auto hit: hitNumbers)
    {
        sum += hit.second;
    }
    cout << sum << endl;
    if(sum < 1)
    {
        
        this->finalResult[this->notConsidered] ++;
        return;
    }
    
    auto krakenShort = this->pruinedTree->getTheMaximumKRAKENhit(hitNumbers);
    
    
    auto UID = this->pruinedTree->getTheUIDFromShort(krakenShort);
    
    auto speciesUID = this->bigTree->getGenusUID(UID);
  
    if(speciesUID == -1)
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

map<LONGS, int>  Tester::getUnieativeHitsWeightedSpecies(YRJObject * yrj , int differences )
{
    map<LONGS, int> unieativeRet;
    
    for(auto kmer : yrj->kmersVector)
    {
        auto hits = this->hits_kmer_with_differences_weighted(kmer, differences);
        
        if(hits.size() == 0){
            continue;
        }
        
        map<LONGS , int> tempHits;
        
        for(auto hit : hits)
        {
            auto speciesUID = getSpeciesLevelUID(hit.first);
            
            if(speciesUID != -1)
            {
                if(tempHits.count(speciesUID))
                    tempHits[speciesUID] = max(tempHits[speciesUID] , hit.second);
                else
                    tempHits[speciesUID] = hit.second;
            }
            else
                cerr << "big problem 44 \n" << hit.first <<endl;
        }
        
        for(auto tmpHit : tempHits)
        {
            if(unieativeRet.count(tmpHit.first))
                unieativeRet[tmpHit.first] += tmpHit.second;
            else
                unieativeRet[tmpHit.first] = tmpHit.second;
        }

    }
    
    return unieativeRet;
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
            auto speciesUID = getSpeciesLevelUID(hit);
            
            if(genusUID == -1 && speciesUID != -1)
            {
                cout << "oooooh " << speciesUID << endl;
                ofstream str("/export1/project/hondius/testingUnieative/newResults/log.txt",std::ofstream::out | std::ofstream::app);
                str << speciesUID << endl;
                continue;
            }
            
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


void Tester::calculate_accurcy_matlab(string file , unordered_map<string , long> finalResult)
{
    ofstream accurFile(file);
    
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
    
    
    speciesAccuracy += getElementInTheResult("species" , finalResult);
    speciesAccuracy += getElementInTheResult("subspecies" , finalResult);
    
    
    
    
    double genusAccuracy = speciesAccuracy ;
    genusAccuracy +=getElementInTheResult ("genus" , finalResult);
    
    
    
    auto family = genusAccuracy +  getElementInTheResult ("family" , finalResult);
    auto order = family +  getElementInTheResult ("order" , finalResult);
    auto classAcc = order +  getElementInTheResult ("class" , finalResult);
    auto phylum = classAcc +  getElementInTheResult ("phylum" , finalResult);
    auto root = phylum +  getElementInTheResult ("root" , finalResult);

    
    accurFile << speciesAccuracy << endl;
    accurFile << genusAccuracy - speciesAccuracy << endl;
    accurFile << family - genusAccuracy  << endl;
    accurFile << order - family   << endl;
    accurFile <<  classAcc - order  << endl;
    accurFile << phylum - classAcc  << endl;
    accurFile << root - phylum  << endl;
    
    
    accurFile.close();
    
    ofstream newFile(file + "garbage");
    
    for(auto ele: finalResult)
        newFile << ele.first << "  " << ele.second << endl;
    
    newFile.close();
    
}



void Tester::calculate_accurcy(string file , unordered_map<string , long> finalResult)
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
    
    
    speciesAccuracy += getElementInTheResult("species", finalResult);
    speciesAccuracy += getElementInTheResult("subspecies" , finalResult);
    
    
    
    
    double genusAccuracy = speciesAccuracy ;
    genusAccuracy +=getElementInTheResult ("genus" , finalResult);
    
    
    double not_counted = getElementInTheResult(notConsidered , finalResult) + getElementInTheResult(notNeeded , finalResult);
    
    
    output << "species:\t" << speciesAccuracy/total << endl;
    output << "genus:\t" << genusAccuracy/total << endl;
    output << "not used\t" << not_counted/(total + not_counted) << endl;
    
    output.close();
    
    
    
    
}



long Tester::getElementInTheResult(string element , unordered_map<string , long> finalResult)
{
    if(finalResult.count(element) )
        return finalResult[element] ;

    return 0;
}




void Tester::testKrakenOutput()
{
    ifstream krakenOutput(path_to_Kraken_test);
    
    string line;
    
    while (getline(krakenOutput, line))
    {
        stringstream converter(line);
        long uid1 , uid2;
        
        converter >> uid1;
        converter >> uid2;
        
        cout << uid1 << "   " << uid2 << endl;
        
        
        auto speciesUID = this->bigTree->getGenusUID(uid2);
        
        if(speciesUID == -1)
        {
            finalResult[this->notNeeded]++;
            continue;
        }

        speciesUID = this->bigTree->getGenusUID(uid1);
        
        if(speciesUID == -1)
        {
            finalResult[this->notNeeded]++;
            continue;
        }
        

        
        auto node1 = bigTree->getNodeFromIndex(bigTree->uid_to_index(uid1));
        cout << "get2\n";
        auto node2 = bigTree->getNodeFromIndex(bigTree->uid_to_index(uid2));
        
        
        
        
        cout << "get3\n";

        
        auto lca = bigTree->get_LCA_between_Two_Nodes(node1, node2);
        cout << "get4\n";

        
        auto  level = this->bigTree->getNextNotNoLevellevel(bigTree->getNodeFromIndex((lca)));
        
        
        cout << "get6\n";

        
        if(this->finalResult.count(level))
            finalResult[level] ++;
        else
            finalResult[level] = 1;


    }
    
}


























void Tester::testKmerLevelLevelMaxMin(YRJObject * yrj)
{
    auto yrjUID = yrj->uid;
    
    
    for(auto kmer : yrj->kmersVector)
    {
        vector<vector<short> > hitsVectors(16);
        auto hits = getHitsandDifferencesKmer(kmer);
        
        for(auto hit : hits)
            hitsVectors[hit.second].emplace_back(hit.first);
        
        // to test this kmer
        
        for (int i= 0  , n = (int)hitsVectors.size() ; i < n  ; ++i)
        {
            if(hitsVectors[i].size() == 0) continue;
            auto pairMaxMin =  getMaxMinLevels(yrj , hitsVectors[i]);
            
            if(max_min_final_results[i].first.count(pairMaxMin.first))
                max_min_final_results[i].first[pairMaxMin.first]++;
            else
                max_min_final_results[i].first[pairMaxMin.first] = 1;
            
            if(max_min_final_results[i].second.count(pairMaxMin.second))
                max_min_final_results[i].second[pairMaxMin.second]++;
            else
                max_min_final_results[i].second[pairMaxMin.second] = 1;
            
            //if(hitsVectors[i].size() > 0) break; // this for purification
        }
    }
}






void Tester::testKmerLevelLevel(YRJObject * yrj)
{
   
    auto yrjUID = yrj->uid;
    
    for(auto kmer : yrj->kmersVector)
    {
        vector<vector<short> > hitsVectors(16);

        auto hits = getHitsandDifferencesKmer(kmer);
        
        for(auto hit : hits)
                hitsVectors[hit.second].emplace_back(hit.first);
    
        
        
        
        // to test this kmer
        
        for (int i= 0  , n = hitsVectors.size() ; i < n  ; ++i)
        {
            
            auto lca = pruinedTree->getGlobalLCA(hitsVectors[i]);
            cout << lca << endl;
            if(lca > -1)
            {
                auto lcaUID = pruinedTree->getTheUIDFromShort(lca);
                cout << lcaUID << endl;
                cout << yrjUID << endl;
                auto level = getLevel(yrjUID, lcaUID);
                
                
                if(finalKmerResult[i].count(level))
                    finalKmerResult[i][level]++;
                else
                    finalKmerResult[i][level] = 1;
                
                //break; // for the pure kmers
            }
        }
    }

    
    

}


map<  short, short > Tester::getHitsandDifferencesKmer(LONG kmer)
{
    map<short, short> hits;
    
    
        for(auto hash : hashes)
        {
            auto tempHits = getNumerOfDifferences(hash, kmer);
            
            for( auto hit : tempHits)
            {
                if(hits.count(hit.first))
                    hits[hit.first] = min(hits[hit.first] , hit.second);
                else
                    hits[hit.first] = hit.second;
            }
            
        }
    
    
    return hits;

}




string Tester::getLevel(long lcaUID , long inputUID)
{
    cout << "1\n";
    auto lcaNode = bigTree->getNodeFromIndex(bigTree->uid_to_index(lcaUID));
    
    cout << "2\n";
    auto inputNode = bigTree->getNodeFromIndex(bigTree->uid_to_index(inputUID));
    
    
    cout << "3\n";
    auto finalLCA =   bigTree->get_LCA_between_Two_Nodes(lcaNode, inputNode);
    
    
    //cout << "4\n";
    
//auto bigIndex = bigTree->uid_to_index(finalLCA);
    
    cout << "5\n";
    auto node = bigTree->getNodeFromIndex(finalLCA);
    
    cout << "6\n";
    return bigTree->getNextNotNoLevellevel(node);
}






void  Tester::writeTestKmerLevels()
{
    
    
    int i = 0;
    for (auto finalMap : finalKmerResult)
    {
        string s;
        if (i < 10)
            s.push_back( (char)(i++ + '0'));
        else
        {
            s.push_back('1');
            s.push_back(i++%10 + '0');
            
        }
        
        calculate_accurcy(this->result + "_summery_" + s + ".out" , finalMap);
        calculate_accurcy_matlab(this->result + "_summery_matlab_" + s + ".out" , finalMap);
        
        
        string ranks[] = {"root" , "phylum" , "class" , "order" , "family" , "genus" , "species" , "subspecies" , "no" , notConsidered , notNeeded};
        
        
        ofstream  *result = new ofstream(this->result + "_summery_" + s + ".out");
        for(auto rank : ranks)
        {
            *result << rank  << "\t" << getElementInTheResult(rank , finalMap) << endl;
        }
        
        result->close();
    }

    
    
    
}












pair<string, string> Tester::getMaxMinLevels(YRJObject * yrj , vector<short> hits)
{
    pair<string, string> max_min_ret;
    max_min_ret.first = "";
    max_min_ret.second = "";
    if(hits.size() == 0)
        return max_min_ret;
    
    
    set<string> results;
    
    for (auto hit : hits)
    {
        results.insert(getLeastCommonLevel(yrj , hit));
    }
    
    
    pair<int, int> max_min_values;
    max_min_values.first  = -1000000;
    max_min_values.second =  1000000;
    
    
    for(auto result : results)
    {
        if(levelsValues[result] > max_min_values.first)
        {
            
            max_min_values.first = levelsValues[result];
            max_min_ret.first = result;
        }
        
        if(levelsValues[result] < max_min_values.second)
        {
            max_min_values.second = levelsValues[result];
            max_min_ret.second = result;
        }
    }
    
    return max_min_ret;
    
}















