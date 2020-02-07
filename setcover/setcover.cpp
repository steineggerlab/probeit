#include "bitset.h"
#include "IITree.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <vector>
#include <algorithm>
#include <limits.h>
struct Probe{
    Probe(int genome, int start, int end, Bitset * set)
     :genome(genome), start(start), end(end), set(set) {}
    int genome;
    int start;
    int end;
    Bitset * set;
};

std::vector<int> setCover(std::vector<Probe> & sets, size_t totalGenomes, int minCovered, int range){
    std::vector<int> coverageCount(totalGenomes);
    std::fill(coverageCount.begin(), coverageCount.end(), 0);
    Bitset covered(totalGenomes);
    std::vector<int> probes;
    std::vector<std::pair<long long, long long>> ranges;
    int max=-1;
    int pos=-1;
    for(size_t i = 0; i < sets.size(); i++){
        int cnt = sets[i].set->Count();
        pos = (cnt > max) ? i : pos;
        max = (cnt > max) ? cnt : max;
    }
    std::cerr << "Total Genomes: " << totalGenomes << std::endl;
    size_t cnt = 0;
    while(covered.Count() < totalGenomes){
        std::cerr << ++cnt <<"\t" << covered.Count() << std::endl;
        probes.push_back(pos);
        for(size_t i = 0; i < totalGenomes; i++){
            if(sets[pos].set->isSet(i)){
                coverageCount[i]++;
                if(coverageCount[i] == minCovered){
                    covered.Set(i);
                }
            }
        }
        int probLen = sets[pos].end - sets[pos].start;
        size_t genomeOffset = (INT_MAX * static_cast<size_t>(sets[pos].genome));

        ranges.emplace_back(genomeOffset + sets[pos].start - (probLen - 1), genomeOffset + sets[pos].end + (probLen - 1));
        IITree<long long, long long> ignore;
        for(size_t i = 0; i < ranges.size(); i++){
            ignore.add(ranges[i].first, ranges[i].second, 0);
        }
        ignore.index();
        //covered.Or(sets[pos]);
        max=-1;
        pos=-1;
        for(size_t i = 0; i < sets.size(); i++){
            int cnt = sets[i].set->AndNotCount(covered);
            size_t genomeOffset = (INT_MAX * static_cast<size_t>(sets[i].genome));
            cnt = (ignore.overlap(genomeOffset + sets[i].start, genomeOffset + sets[i].end)) ? 0 :  cnt;
            pos = (cnt > max) ? i : pos;
            max = (cnt > max) ? cnt : max;
        }
    }
    return probes;
}

std::vector<Probe> readInSet(std::string filePath, int probeLen){
    std::set<int> genomeIds;
    std::fstream infile;
    infile.open(filePath, std::fstream::in);
    std::vector<Probe> genomicSets;
    std::string pipes = "|";
    std::string line;
    while (std::getline(infile, line)) {
        std::string kmerInGenome = line.substr(line.find(";") + 1, line.length());
        size_t pos = 0;
        std::string genomeAndPos;
        while ((pos = kmerInGenome.find(pipes)) != std::string::npos ) {
            genomeAndPos = kmerInGenome.substr(0, pos);
            size_t columnPos = genomeAndPos.find(",");
            int currGenome = stoi(genomeAndPos.substr(0, columnPos));
            genomeIds.insert(currGenome);
            kmerInGenome.erase(0, pos + pipes.length());
        }
    }
    infile.clear();
    infile.seekg(0, std::ios::beg);
    
    size_t totalGenomes = genomeIds.size();
    std::cerr << "Parse" << std::endl;
    size_t simKmerCnt = 0;
    while (std::getline(infile, line)) {
        Bitset * lineBitSet = new Bitset(totalGenomes);
        //0,7;0,7|93,159|1656,227
        std::string currGenomeAndPos = line.substr(0, line.find(";"));
        size_t columnPos = currGenomeAndPos.find(",");
        int genome = stoi(currGenomeAndPos.substr(0, columnPos));
        int genomePos = stoi(currGenomeAndPos.substr(columnPos+1, currGenomeAndPos.length()));
        std::string kmerInGenome = line.substr(line.find(";") + 1, line.length());
        //std::cout <<kmerInGenome << std::endl;
        size_t pos = 0;
        std::string genomeAndPos;
        while ((pos = kmerInGenome.find(pipes)) != std::string::npos ) {
            genomeAndPos = kmerInGenome.substr(0, pos);
            size_t columnPos = genomeAndPos.find(",");
            int currGenome = stoi(genomeAndPos.substr(0, columnPos));
            lineBitSet->Set(currGenome);
            // iterate pipe by pipe
            kmerInGenome.erase(0, pos + pipes.length());
        }
        
        genomeAndPos = kmerInGenome.substr(0, pos);
        lineBitSet->Set(stoi(genomeAndPos.substr(0, genomeAndPos.find(","))));
        genomicSets.emplace_back(genome, genomePos, genomePos+probeLen, lineBitSet);
    }
    infile.close();
    return genomicSets;
}



int main(int argc, char ** argv){
    int minCovered = 1;
    int range = 1;
    if(argc >= 3){
        minCovered = atoi(argv[2]);
    }
    if(argc == 4){
        range = atoi(argv[3]);
    }
    std::vector<Probe> probeSet = readInSet(std::string(argv[1]), 40);
    size_t totalGenomes = probeSet[0].set->Size();

    std::cerr << "SetCover Mincover=" << minCovered << std::endl;
    std::vector<int> cover = setCover(probeSet, totalGenomes, minCovered, range);
    std::sort(cover.begin(), cover.end(), std::greater<int>());
    for(size_t i = 0; i< cover.size(); i++){
        std::cerr << cover[i] << std::endl;
    }

    int cnt = 0;
    int coverLine = cover.back();
    cover.pop_back();
    std::fstream infile;
    infile.open(argv[1], std::fstream::in);
    std::string line;
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
    
    
    for(size_t set = 0; set < probeSet.size(); set++){
        delete probeSet[set].set;
    }
    
    return 0;
}

