#include "bitset.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <vector>
#include <algorithm>

std::vector<int> setCover(std::vector<Bitset> & sets, size_t totalGenomes, int minCovered, int range){
    std::vector<int> coverageCount(totalGenomes);
    std::fill(coverageCount.begin(), coverageCount.end(), 0);
    Bitset covered(totalGenomes);
    std::vector<int> probes;
    int max=-1;
    int pos=-1;
    for(size_t i = 0; i < sets.size(); i++){
        int cnt = sets[i].Count();
        pos = (cnt > max) ? i : pos;
        max = (cnt > max) ? cnt : max;
    }
    std::cerr << "Total Genomes: " << totalGenomes << std::endl;
    size_t cnt = 0;
    while(covered.Count() < totalGenomes){
        std::cerr << ++cnt <<"\t" << covered.Count() << std::endl;
        probes.push_back(pos);
        for(size_t i = 0; i < totalGenomes; i++){
            if(sets[pos].isSet(i)){
                coverageCount[i]++;
                if(coverageCount[i] == minCovered){
                    covered.Set(i);
                }
            }
        }
        for(size_t i = pos; i <= pos+range && i < sets.size(); i++)
            sets[i].Clear();
        for(size_t i = std::max(0, pos-range); i < pos; i++)
            sets[i].Clear();
        //covered.Or(sets[pos]);
        max=-1;
        pos=-1;
        for(size_t i = 0; i < sets.size(); i++){
            int cnt = sets[i].AndNotCount(covered);
            pos = (cnt > max) ? i : pos;
            max = (cnt > max) ? cnt : max;
        }
    }
    return probes;
}


int main(int argc, char ** argv){
    std::set<int> genomeIds;
    std::fstream infile;
    infile.open( std::string(argv[1]), std::fstream::in);
    int minCovered = 1;
    if(argc == 3){
        minCovered = atoi(argv[2]);
    }
    // skip first line
    std::string line;
    std::getline(infile, line);
    while (std::getline(infile, line)) {
        std::string genomeId = line.substr(0, line.find(","));
        genomeIds.insert(stoi(genomeId));
    }
    infile.clear();
    infile.seekg(0, std::ios::beg);
    size_t totalGenomes = genomeIds.size();
    std::vector<Bitset> genomicSets;
    std::cerr << "Parse" << std::endl;
    // skip first line
    std::getline(infile, line);
    size_t simKmerCnt = 0;
    while (std::getline(infile, line)) {
        Bitset lineBitSet(totalGenomes);
        std::string pipes = "|";
        //0,7;0,7|93,159|1656,227
        std::string kmerInGenome = line.substr(line.find(";") + 1, line.length());
        std::string currGenomeAndPos = line.substr(0, line.find(";"));
        size_t columnPos = currGenomeAndPos.find(",");
        int genome = stoi(currGenomeAndPos.substr(0, columnPos));
        int genomePos = stoi(currGenomeAndPos.substr(columnPos+1, currGenomeAndPos.length()));

        //std::cout <<kmerInGenome << std::endl;
        size_t pos = 0;
        std::string genomeAndPos;
        bool sameKmer = false;
        while ((pos = kmerInGenome.find(pipes)) != std::string::npos && sameKmer == false) {
            genomeAndPos = kmerInGenome.substr(0, pos);
            columnPos = genomeAndPos.find(",");
            int currGenome = stoi(genomeAndPos.substr(0, columnPos));
            lineBitSet.Set(currGenome);
            int currGenomePos = stoi(genomeAndPos.substr(columnPos+1, genomeAndPos.length()));
            if(currGenome < genome || (currGenome==genome && currGenomePos < genomePos )){
                sameKmer = true;
                simKmerCnt++;
                lineBitSet.Clear();
            }
            kmerInGenome.erase(0, pos + pipes.length());
        }
        
        genomeAndPos = kmerInGenome.substr(0, pos);
        lineBitSet.Set(stoi(genomeAndPos.substr(0, genomeAndPos.find(","))));
        genomicSets.push_back(lineBitSet);
    }
    std::cerr << "Ignore " << simKmerCnt << " kmers" << std::endl;
    std::cerr << "SetCover Mincover=" << minCovered << std::endl;
    std::vector<int> cover = setCover(genomicSets, totalGenomes, minCovered, 1);
    std::sort(cover.begin(), cover.end(), std::greater<int>());
    for(size_t i = 0; i< cover.size(); i++){
        std::cerr << cover[i] << std::endl;
    }
    infile.clear();
    infile.seekg(0, std::ios::beg);
    int cnt = 0;
    int coverLine = cover.back();
    cover.pop_back();
    while (std::getline(infile, line)) {
        if( cnt == coverLine){
           std::cout << line << std::endl;
           if(cover.size()==0){
              break;
           }
           coverLine = cover.back();
           cover.pop_back();
        }
        cnt++; 
    }
    return 0;
}

