#include <string.h>
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
    Probe(int genome, int start, int end, std::vector<std::pair<int,int>> * genomeAndPosition)
     :genome(genome), start(start), end(end), genomeAndPosition(genomeAndPosition) {}
    int genome;
    int start;
    int end;
    std::vector<std::pair<int,int>> * genomeAndPosition;
};

std::vector<int> setCover(std::vector<Probe> & sets, size_t totalGenomes, int minCovered, int range){
    std::vector<int> coverageCount(totalGenomes);
    std::fill(coverageCount.begin(), coverageCount.end(), 0);
    char * covered = new char [totalGenomes];
    memset(covered, 0, sizeof(char) * totalGenomes);
    int coveredCount = 0;
    std::vector<int> probes;
    std::vector<std::pair<long long, long long>> ranges;
    int maxProbeCover=-1;
    int maxProbeIdx=-1;
    for(size_t i = 0; i < sets.size(); i++){
        int cnt = sets[i].genomeAndPosition->size();
        maxProbeIdx = (cnt > maxProbeCover) ? i : maxProbeIdx;
        maxProbeCover = (cnt > maxProbeCover) ? cnt : maxProbeCover;
    }
    std::cerr << "Total Genomes: " << totalGenomes << std::endl;
    size_t cnt = 0;
    while(coveredCount < totalGenomes && maxProbeCover > 0){
        probes.push_back(maxProbeIdx);
        for(size_t setElmIdx = 0; setElmIdx < sets[maxProbeIdx].genomeAndPosition->size(); setElmIdx++){
            int elmIdx = sets[maxProbeIdx].genomeAndPosition->at(setElmIdx).first;
            coverageCount[elmIdx]++;
            if(coverageCount[elmIdx] == minCovered){
                covered[elmIdx] = 1;
                coveredCount++;
            }
        }
        std::cerr << ++cnt <<"\t" << coveredCount << std::endl;

        int probLen = sets[maxProbeIdx].end - sets[maxProbeIdx].start;
        IITree<long long, long long> ignore;
        for(size_t setElmIdx = 0; setElmIdx < sets[maxProbeIdx].genomeAndPosition->size(); setElmIdx++){
            int genome   = sets[maxProbeIdx].genomeAndPosition->at(setElmIdx).first;
            int startPos = sets[maxProbeIdx].genomeAndPosition->at(setElmIdx).second;
            size_t genomeOffset = (INT_MAX * static_cast<size_t>(genome));
            ranges.emplace_back(genomeOffset + startPos - (probLen - 1), genomeOffset + startPos + probLen + 1);
        }
        for(size_t i = 0; i < ranges.size(); i++){
            ignore.add(ranges[i].first, ranges[i].second, 0);
        }
        ignore.index();
        //covered.Or(sets[pos]);
        maxProbeCover=-1;
        maxProbeIdx=-1;
        for(size_t i = 0; i < sets.size(); i++) {
            int setCoverElement = 0;
            for(size_t setElmIdx = 0; setElmIdx < sets[i].genomeAndPosition->size(); setElmIdx++){
                int genome = sets[i].genomeAndPosition->at(setElmIdx).first;
                setCoverElement += (covered[genome]);
            }
            
            int cnt = sets[i].genomeAndPosition->size() - setCoverElement;
            bool hasOverlap = false;
            for(size_t setElmIdx = 0; setElmIdx < sets[i].genomeAndPosition->size(); setElmIdx++) {
                int genome = sets[i].genomeAndPosition->at(setElmIdx).first;
                int startPos = sets[i].genomeAndPosition->at(setElmIdx).second;
                size_t genomeOffset = (INT_MAX * static_cast<size_t>(genome));
                hasOverlap |= ignore.overlap(genomeOffset + startPos - (probLen - 1), genomeOffset + startPos + probLen + 1);
            }
            cnt = (hasOverlap) ? 0 :  cnt;
            maxProbeIdx = (cnt > maxProbeCover) ? i : maxProbeIdx;
            maxProbeCover = (cnt > maxProbeCover) ? cnt : maxProbeCover;
        }
    }
    return probes;
}

std::vector<Probe> readInSet(std::string filePath, size_t & totalGenomes, int probeLen){
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
    
    totalGenomes = genomeIds.size()+1;
    std::cerr << "Parse" << std::endl;
    while (std::getline(infile, line)) {
        //0,7;0,7|93,159|1656,227
        std::string currGenomeAndPos = line.substr(0, line.find(";"));
        size_t columnPos = currGenomeAndPos.find(",");
        int genome = stoi(currGenomeAndPos.substr(0, columnPos));
        int genomePos = stoi(currGenomeAndPos.substr(columnPos+1, currGenomeAndPos.length()));
        std::string kmerInGenome = line.substr(line.find(";") + 1, line.length());
        //std::cout <<kmerInGenome << std::endl;
        size_t pos = 0;
        std::string genomeAndPos;
        std::vector<std::pair<int, int>> * genomeAndPosition = new std::vector<std::pair<int, int>>();
        while ((pos = kmerInGenome.find(pipes)) != std::string::npos ) {
            genomeAndPos = kmerInGenome.substr(0, pos);
            size_t columnPos = genomeAndPos.find(",");
            int currGenome = stoi(genomeAndPos.substr(0, columnPos));
            int currGenomePos = stoi(genomeAndPos.substr(columnPos+1, genomeAndPos.length()));
            genomeAndPosition->emplace_back(currGenome, currGenomePos);
            // iterate pipe by pipe
            kmerInGenome.erase(0, pos + pipes.length());
        }
        
        genomeAndPos = kmerInGenome.substr(0, pos);
        genomicSets.emplace_back(genome, genomePos, genomePos+probeLen, genomeAndPosition);
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
    size_t totalGenomes;
    std::vector<Probe> probeSet = readInSet(std::string(argv[1]), totalGenomes, 40);

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
        if(cnt == coverLine){
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
        delete probeSet[set].genomeAndPosition;
    }
    
    return 0;
}

