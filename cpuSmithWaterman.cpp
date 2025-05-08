#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cctype>
#include <cstdio>
#include <algorithm>
#include <chrono>  // For timing
#include <iomanip>  // For std::setprecision

// Helper function to extract base name without directory or extension
std::string extractBaseName(const std::string& filepath) {
    // Find the last slash or backslash
    size_t lastSlash = filepath.find_last_of("/\\");
    size_t start = (lastSlash == std::string::npos) ? 0 : lastSlash + 1;
    
    // Find the last dot (extension)
    size_t lastDot = filepath.find_last_of('.');
    size_t end = (lastDot == std::string::npos || lastDot <= start) ? filepath.length() : lastDot;
    
    // Extract the filename without extension
    return filepath.substr(start, end - start);
}

// Strip any prefix like "BB11001_" from sequence names
std::string stripPrefix(const std::string& name) {
    size_t pos = name.find('_');
    if(pos != std::string::npos && pos > 0) {
        return name.substr(pos + 1);
    }
    return name;
}

// Read a sequence from a FASTA file
bool readFastaFile(const std::string& filename, std::string& name, std::string& seq) {
    std::ifstream fin(filename);
    if(!fin.is_open()) {
        return false;
    }
    
    std::string line;
    name = "";
    seq = "";
    
    // Read header line
    if(std::getline(fin, line)) {
        if(line.size() > 0 && line[0] == '>') {
            // Extract name (up to first whitespace or end of line after '>')
            size_t pos = 1;
            while(pos < line.size() && !isspace(static_cast<unsigned char>(line[pos]))) {
                name.push_back(line[pos]);
                pos++;
            }
        }
        
        // Read sequence lines
        while(std::getline(fin, line)) {
            if(line.size() > 0 && line[0] == '>') break; // stop if another header
            for(char c : line) {
                if(!isspace(static_cast<unsigned char>(c))) {
                    seq.push_back(c);
                }
            }
        }
    }
    
    fin.close();
    return true;
}

// Convert a sequence to uppercase
void toUpperCase(std::string& seq) {
    for(char &c : seq) {
        c = std::toupper(static_cast<unsigned char>(c));
    }
}

// Perform Smith-Waterman alignment
void smithWaterman(const std::string& seq1, const std::string& seq2, 
                  int matchScore, int mismatchScore, int gapScore,
                  std::string& align1, std::string& align2, int& maxScore) {
    int len1 = seq1.length();
    int len2 = seq2.length();
    
    // Allocate score and direction matrices
    std::vector<int> score((len1+1) * (len2+1), 0);
    std::vector<unsigned char> dir((len1+1) * (len2+1), 0);
    
    // Fill the matrices
    for(int i = 1; i <= len1; ++i) {
        for(int j = 1; j <= len2; ++j) {
            // Compute scores for match/mismatch and gap options
            int up   = score[(i-1) * (len2+1) + j] + gapScore;
            int left = score[i * (len2+1) + (j-1)] + gapScore;
            int diagScore = score[(i-1) * (len2+1) + (j-1)] + 
                           ((seq1[i-1] == seq2[j-1]) ? matchScore : mismatchScore);
            
            // Choose the maximum, compare with 0 for local alignment
            int localMaxScore = 0;
            unsigned char direction = 0;
            if(diagScore > localMaxScore) {
                localMaxScore = diagScore;
                direction = 1; // 1 = diagonal
            }
            if(up > localMaxScore) {
                localMaxScore = up;
                direction = 2; // 2 = up (gap in seq2)
            }
            if(left > localMaxScore) {
                localMaxScore = left;
                direction = 3; // 3 = left (gap in seq1)
            }
            
            // Store score and direction
            score[i * (len2+1) + j] = localMaxScore;
            dir[i * (len2+1) + j] = direction;
        }
    }
    
    // Find the cell with maximum score
    maxScore = 0;
    int max_i = 0, max_j = 0;
    for(int i = 1; i <= len1; ++i) {
        for(int j = 1; j <= len2; ++j) {
            int val = score[i * (len2+1) + j];
            if(val > maxScore) {
                maxScore = val;
                max_i = i;
                max_j = j;
            }
        }
    }
    
    // Traceback from (max_i, max_j) until score becomes 0
    align1 = "";
    align2 = "";
    int ti = max_i;
    int tj = max_j;
    while(ti > 0 && tj > 0) {
        unsigned char d = dir[ti * (len2+1) + tj];
        if(d == 0) {
            break; // alignment stop
        }
        if(d == 1) { // diagonal
            align1.push_back(seq1[ti-1]);
            align2.push_back(seq2[tj-1]);
            ti -= 1;
            tj -= 1;
        } else if(d == 2) { // came from up (gap in seq2)
            align1.push_back(seq1[ti-1]);
            align2.push_back('-');  // use '-' for gap during traceback
            ti -= 1;
        } else if(d == 3) { // came from left (gap in seq1)
            align1.push_back('-');
            align2.push_back(seq2[tj-1]);
            tj -= 1;
        } else {
            // Should not happen for Smith-Waterman (d is 0-3)
            break;
        }
        if(score[ti * (len2+1) + tj] == 0) {
            // Stop when we hit a cell with 0 (beginning of local alignment)
            break;
        }
    }
    
    // Reverse the aligned strings as we collected them backward
    std::reverse(align1.begin(), align1.end());
    std::reverse(align2.begin(), align2.end());
    
    // Replace '-' gaps with '.' for MSF output format
    for(char &c : align1) {
        if(c == '-') c = '.';
    }
    for(char &c : align2) {
        if(c == '-') c = '.';
    }
}

// Compute GCG checksum for a sequence
int gcgChecksum(const std::string &s) {
    long check = 0;
    for(size_t i = 0; i < s.size(); ++i) {
        // use (i % 57) + 1 multiplier
        int c = std::toupper(static_cast<unsigned char>(s[i]));
        check += ((int)((i % 57) + 1) * c);
    }
    return (int)(check % 10000);
}

// Determine sequence type (DNA or protein)
char determineSequenceType(const std::string& align1, const std::string& align2) {
    std::string allAligned = align1 + align2;
    bool maybeDNA = true;
    for(char c : allAligned) {
        // skip gaps
        if(c == '.') continue;
        char u = std::toupper(static_cast<unsigned char>(c));
        if(u != 'A' && u != 'C' && u != 'G' && u != 'T' && u != 'U' && u != 'N') {
            maybeDNA = false;
            break;
        }
    }
    return maybeDNA ? 'N' : 'P';
}

// Print alignment in MSF format
void printMSFAlignment(const std::string& name1, const std::string& name2,
                     const std::string& align1, const std::string& align2,
                     int maxScore) {
    int alignLen = align1.size();
    
    // Compute checksums
    int check1 = gcgChecksum(align1);
    int check2 = gcgChecksum(align2);
    int globalCheck = ((long)check1 + (long)check2) % 10000;
    
    // Determine sequence type
    char typeChar = determineSequenceType(align1, align2);
    
    std::cout << "Alignment score: " << maxScore << "\n\n";
    
    // Output alignment in MSF (PileUp) format
    std::cout << "PileUp\n\n";
    std::printf("   MSF:   %d  Type: %c    Check:  %4d   ..\n\n", alignLen, typeChar, globalCheck);
    std::printf(" Name: %s oo  Len:   %d  Check:  %4d  Weight:  10.0\n", name1.c_str(), alignLen, check1);
    std::printf(" Name: %s oo  Len:   %d  Check:  %4d  Weight:  10.0\n\n", name2.c_str(), alignLen, check2);
    std::cout << "//\n\n";
    
    // Print aligned sequences in blocks of 50 columns
    int colsPerLine = 50;
    for(int start = 0; start < alignLen; start += colsPerLine) {
        int end = (start + colsPerLine < alignLen) ? (start + colsPerLine) : alignLen;
        // Sequence 1 line
        std::printf("%-12s", name1.c_str());  // name left padded to 12 characters
        // Print sequence with a space every 10 residues
        int count = 0;
        for(int k = start; k < end; ++k) {
            std::cout << align1[k];
            count++;
            if(count % 10 == 0 && k < end - 1) {
                std::cout << ' ';
            }
        }
        std::cout << "\n";
        // Sequence 2 line
        std::printf("%-12s", name2.c_str());
        count = 0;
        for(int k = start; k < end; ++k) {
            std::cout << align2[k];
            count++;
            if(count % 10 == 0 && k < end - 1) {
                std::cout << ' ';
            }
        }
        std::cout << "\n\n";
    }
}

int main(int argc, char **argv) {
    // Start timing the execution
    auto startTime = std::chrono::high_resolution_clock::now();
    
    if(argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <seq1.fasta> <seq2.fasta>\n";
        return 1;
    }
    std::string file1 = argv[1];
    std::string file2 = argv[2];
    
    // Read sequences from FASTA files
    std::string name1, name2;
    std::string seq1, seq2;
    
    if(!readFastaFile(file1, name1, seq1) || !readFastaFile(file2, name2, seq2)) {
        std::cerr << "Error: unable to open or parse input FASTA file(s).\n";
        return 1;
    }
    
    // If names are not provided in FASTA, use filenames instead
    if(name1.empty()) {
        name1 = extractBaseName(file1);
    }
    if(name2.empty()) {
        name2 = extractBaseName(file2);
    }
    
    // Strip any prefixes from sequence names
    name1 = stripPrefix(name1);
    name2 = stripPrefix(name2);
    
    // Convert sequences to uppercase
    toUpperCase(seq1);
    toUpperCase(seq2);
    
    if(seq1.empty() || seq2.empty()) {
        std::cerr << "Error: one of the sequences is empty.\n";
        return 1;
    }
    
    // Scoring scheme
    int matchScore = 2;
    int mismatchScore = -1;
    int gapScore = -1;
    
    // Perform Smith-Waterman alignment
    std::string align1, align2;
    int maxScore;
    
    smithWaterman(seq1, seq2, matchScore, mismatchScore, gapScore, 
                  align1, align2, maxScore);
    
    // Calculate and output execution time with microsecond precision
    auto endTime = std::chrono::high_resolution_clock::now();
    
    // Calculate durations in different units for more precise reporting
    auto durationMicro = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
    auto durationNano = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count();
    
    // Output timing in the most appropriate unit
    if (durationMicro < 10000) {  // Less than 10ms, show in microseconds
        std::cerr << "CPU Execution time: " << durationMicro << " μs (" << durationNano << " ns)" << std::endl;
    } else {
        // For longer runtimes, show in milliseconds with microsecond precision
        double durationMs = static_cast<double>(durationMicro) / 1000.0;
        std::cerr << "CPU Execution time: " << std::fixed << std::setprecision(3) << durationMs << " ms" << std::endl;
    }
    
    // Print alignment in MSF format
    printMSFAlignment(name1, name2, align1, align2, maxScore);
    
    return 0;
}