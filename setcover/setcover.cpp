#include <string.h>
#include "kseq.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <vector>
#include <algorithm>
#include <limits.h>
#include <unistd.h>

KSEQ_INIT(int, read)

struct Probe{
    Probe(size_t line, int genome, int start, int end, std::vector<std::pair<int,int>> * genomeAndPosition)
     :line(line), genome(genome), start(start), end(end), genomeAndPosition(genomeAndPosition) {}
    size_t line;
    int genome;
    int start;
    int end;
    std::vector<std::pair<int,int>> * genomeAndPosition;

    static bool comparByLine (const Probe &first, const Probe &second) {
        return (first.line > second.line);
    }
};

/*
   * Adapted from levenshtein.js (https://gist.github.com/andrei-m/982927)
   * Changed to hold only one row of the dynamic programing matrix
   * Copyright (c) 2011 Andrei Mackenzie
   * Martin Steinegger: Changed to local alignment version
   * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
   * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
   * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
static int localLevenshteinDistance(const std::string &s1, const std::string &s2) {
    int m = s1.length();
    int n = s2.length();

    if (m == 0)
        return n;
    if (n == 0)
        return m;

    // we don't need the full matrix just for the distance
    // we only need space for one row + the last value;
    int *currRow = new int[m + 1];
    for (int i = 0; i <= m; ++i) {
        currRow[i] = 0;
    }
    int maxVal=0;
    int val;
    // fill in the rest of the matrix
    for (int i = 1; i <= n; i++) {
        int prev = 0;
        for (int j = 1; j <= m; j++) {
            int subScore = (s2[i - 1] == s1[j - 1]) ? currRow[j - 1] + 1 : currRow[j - 1] - 2;
            val = std::max(0,
                           std::max(subScore, // substitution
                                    std::max(prev - 2,       // insertion
                                             currRow[j] - 2)));   // deletion


            if(val > maxVal){
                maxVal = val;
            }
            currRow[j - 1] = prev;
            prev = val;
        }
//            std::swap(currRow, topRow);
    }

    // the last element of the matrix is the distance
    delete[] currRow;

    return maxVal;
}


bool checkAlignment(std::vector<std::string> & seqs, Probe &probe, std::vector<Probe> &pickedProbes, int distanceThreshhold){
    int probeLen = probe.end - probe.start;
    std::string queryProbe(seqs[probe.genome].begin()+probe.start, seqs[probe.genome].begin()+probe.end);
    for(size_t i = 0; i < pickedProbes.size(); i++){
        std::string targetProbe(seqs[pickedProbes[i].genome].begin()+pickedProbes[i].start,
                                seqs[pickedProbes[i].genome].begin()+pickedProbes[i].end);
        int distance = localLevenshteinDistance(queryProbe, targetProbe);
        if(distance > distanceThreshhold){
            return false;
        }
    }

    return true;
}

std::vector<std::string> readFasta(std::string &fastaFile) {
    kseq_t *seq;
    std::vector<std::string> seqs;
    FILE * fp = fopen(fastaFile.c_str(), "r"); // STEP 2: open the file handler
    seq = kseq_init(fileno(fp)); // STEP 3: initialize seq
    while ( kseq_read(seq) >= 0) { // STEP 4: read sequence
        seqs.push_back(seq->seq.s);
    }
    kseq_destroy(seq);
    fclose(fp);
    return seqs;
}

std::vector<Probe> setCover(std::vector<Probe> & sets, std::vector<std::string> & seqs,
                          size_t totalGenomes, int minCovered, float percentCoverd, int distanceThreshhold){
    std::vector<int> coverageCount(totalGenomes);
    std::fill(coverageCount.begin(), coverageCount.end(), 0);
    std::vector<bool> * isSet = new std::vector<bool>[seqs.size()];
    for(size_t i = 0; i < seqs.size(); i++){
        isSet[i].reserve(seqs[i].size());
        std::fill(isSet[i].begin(), isSet[i].begin()+seqs[i].size(), 0);
    }
    char * covered = new char [totalGenomes];
    memset(covered, 0, sizeof(char) * totalGenomes);
    int coveredCount = 0;
    std::vector<Probe> probes;
    std::vector<std::pair<long long, long long>> ranges;
    int maxProbeCover=-1;
    int maxProbeIdx=-1;
    for(size_t i = 0; i < sets.size(); i++){
        int cnt = sets[i].genomeAndPosition->size();
        maxProbeIdx   = (cnt > maxProbeCover) ? i : maxProbeIdx;
        maxProbeCover = (cnt > maxProbeCover) ? cnt : maxProbeCover;
    }
    std::cerr << "Total Genomes: " << totalGenomes << std::endl;
    size_t cnt = 0;
    while(coveredCount < totalGenomes && maxProbeCover > 0){
        if(checkAlignment(seqs, sets[maxProbeIdx], probes, distanceThreshhold)){
            probes.push_back(sets[maxProbeIdx]);
            for(size_t setElmIdx = 0; setElmIdx < sets[maxProbeIdx].genomeAndPosition->size(); setElmIdx++){
                int elmIdx = sets[maxProbeIdx].genomeAndPosition->at(setElmIdx).first;
                coverageCount[elmIdx]++;
                if(coverageCount[elmIdx] == minCovered){
                    covered[elmIdx] = 1;
                    coveredCount++;
                }
            }
            std::cerr << ++cnt <<"\t" << coveredCount << std::endl;
        }
        if(static_cast<float>(coveredCount)/static_cast<float>(totalGenomes) >= percentCoverd ){
            break;
        }

        int probLen = sets[maxProbeIdx].end - sets[maxProbeIdx].start;
        for(size_t setElmIdx = 0; setElmIdx < sets[maxProbeIdx].genomeAndPosition->size(); setElmIdx++){
            int genome   = sets[maxProbeIdx].genomeAndPosition->at(setElmIdx).first;
            int startPos = sets[maxProbeIdx].genomeAndPosition->at(setElmIdx).second;
            for(size_t pos = startPos; pos < startPos+probLen; pos++){
                isSet[genome][pos] = 1;
            }
        }
        //covered.Or(sets[pos]);
        maxProbeCover=0;
        maxProbeIdx=-1;
        for(size_t i = 0; i < sets.size(); i++) {
            int setCoverElement = 0;
            // cnt > sets[i].genomeAndPosition->size() <- should speed up
            // if already a cnt was found that is > the possible cnt
            // than we do not need to evaluate
            if(maxProbeCover >= sets[i].genomeAndPosition->size()){
                continue;
            }
            bool hasOverlap = false;
            for(size_t setElmIdx = 0; setElmIdx < sets[i].genomeAndPosition->size() && hasOverlap == false; setElmIdx++) {
                int genome = sets[i].genomeAndPosition->at(setElmIdx).first;
                int startPos = sets[i].genomeAndPosition->at(setElmIdx).second;
                hasOverlap |= isSet[genome][startPos] | isSet[genome][startPos+probLen-1];
            }
            if(hasOverlap){
                continue;
            }

            for(size_t setElmIdx = 0; setElmIdx < sets[i].genomeAndPosition->size(); setElmIdx++){
                int genome = sets[i].genomeAndPosition->at(setElmIdx).first;
                setCoverElement += (covered[genome]);
            }
            
            int setCnt = sets[i].genomeAndPosition->size() - setCoverElement;
            maxProbeIdx = (setCnt > maxProbeCover) ? i : maxProbeIdx;
            maxProbeCover = (setCnt > maxProbeCover) ? setCnt : maxProbeCover;
        }
    }
    delete [] isSet;
    return probes;
}

std::vector<Probe> readInSet(std::string &filePath, size_t & totalGenomes, int probeLen){
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
        size_t columnPos = kmerInGenome.find(",");
        int currGenome = stoi(kmerInGenome.substr(0, columnPos));
        genomeIds.insert(currGenome);
    }
    infile.clear();
    infile.seekg(0, std::ios::beg);
    std::set<int>::iterator it;
    for (it = genomeIds.begin(); it != genomeIds.end(); ++it)
    {
        totalGenomes = std::max(totalGenomes, static_cast<size_t>(*it)); // Note the "*" here
    }
    totalGenomes++;
    size_t lineCnt = 0;
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
        columnPos = kmerInGenome.find(",");
        int currGenome = stoi(kmerInGenome.substr(0, columnPos));
        int currGenomePos = stoi(kmerInGenome.substr(columnPos+1, kmerInGenome.length()));
        genomeAndPosition->emplace_back(currGenome, currGenomePos);

        genomeAndPos = kmerInGenome.substr(0, pos);
        genomicSets.emplace_back(lineCnt, genome, genomePos, genomePos+probeLen, genomeAndPosition);
        lineCnt++;
    }
    infile.close();
    return genomicSets;
}

int main(int argc, char ** argv){
    int minCovered = 1;
    int range = 1;
    float percentCoverd = 1.0f;
    int distanceThreshhold = UINT_MAX;
    int randIterations = 1;
    if(argc == 1){
        std::cout << "Usage: setcover [options] <genmapTSV> <fasta>\n";
        std::cout << "Options:\n";
        std::cout << "  -c INT\teach sequence should be N times covered by a probe\n";
        std::cout << "  -l INT\tlength of probe. Set if probes should not overlap else set 1\n";
        std::cout << "  -p FLOAT\tearly stop if X% of the sequences are '-c N' times covered\n";
        std::cout << "  -d INT\taccept probes only if the levenshtein distance to previously picked probes is < than -d\n";
        std::cout << "  -i INT\tminimize probes randomly N iterations \n";
        exit(0);
    }
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-c") && (i < argc - 1)) {
            minCovered = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-l") && (i < argc - 1)) {
            range = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-p") && (i < argc - 1)) {
            percentCoverd = atof(argv[++i]);
        } else if (!strcmp(argv[i], "-d") && (i < argc - 1)) {
            distanceThreshhold = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-i") && (i < argc - 1)) {
            randIterations = atoi(argv[++i]);
        }
    }

    std::string fastaFile = std::string(argv[argc-1]);
    std::string genmapTSV = std::string(argv[argc-2]);
    std::vector<std::string> seqs = readFasta(fastaFile);

    size_t totalGenomes = 0;
    std::vector<Probe> probeSet = readInSet(genmapTSV, totalGenomes, range);
    
    std::cerr << "SetCover Mincover=" << minCovered << std::endl;
    int coverProbCnt = INT_MAX;
    std::vector<Probe> result;
    for(size_t i = 0; i < randIterations; i++){
        // std::random_shuffle ( probeSet.begin(), probeSet.end() );
        std::vector<Probe> cover = setCover(probeSet, seqs, totalGenomes, minCovered, percentCoverd, distanceThreshhold);
        if(cover.size() < coverProbCnt){
            coverProbCnt = cover.size();
            result = cover;
        }
    }
    std::sort(result.begin(), result.end(), Probe::comparByLine);

    int cnt = 0;
    int coverLine = result.back().line;
    result.pop_back();
    std::fstream infile;
    infile.open(genmapTSV, std::fstream::in);
    std::string line;
    while (std::getline(infile, line)) {
        if(cnt == coverLine){
           std::cout << line << std::endl;
           if(result.size()==0){
              break;
           }
           coverLine = result.back().line;
           result.pop_back();
        }
        cnt++; 
    }
    
    for(size_t set = 0; set < probeSet.size(); set++){
        delete probeSet[set].genomeAndPosition;
    }
    
    return 0;
}


