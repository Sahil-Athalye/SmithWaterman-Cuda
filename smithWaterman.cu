#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cctype>
#include <cstdio>
#include <algorithm>

// CUDA kernel to compute one anti-diagonal of the Smith-Waterman DP and direction matrices
__global__ void sw_kernel(const char *seq1, const char *seq2, int len1, int len2, 
                          int diag, int start_i, int end_i, 
                          int *score, unsigned char *dir, int matchScore, int mismatchScore, int gapScore) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i = start_i + idx;
    if(i > end_i) return;
    int j = diag - i;
    // Compute scores for match/mismatch and gap options
    int up   = score[(i-1) * (len2+1) + j] + gapScore;
    int left = score[i * (len2+1) + (j-1)] + gapScore;
    int diagScore = score[(i-1) * (len2+1) + (j-1)] + ((seq1[i-1] == seq2[j-1]) ? matchScore : mismatchScore);
    // Choose the maximum, compare with 0 for local alignment
    int maxScore = 0;
    unsigned char direction = 0;
    if(diagScore > maxScore) {
        maxScore = diagScore;
        direction = 1; // 1 = diagonal
    }
    if(up > maxScore) {
        maxScore = up;
        direction = 2; // 2 = up (gap in seq2)
    }
    if(left > maxScore) {
        maxScore = left;
        direction = 3; // 3 = left (gap in seq1)
    }
    // Write back score and direction
    score[i * (len2+1) + j] = maxScore;
    dir[i * (len2+1) + j] = direction;
}

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

int main(int argc, char **argv) {
    if(argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <seq1.fasta> <seq2.fasta>\n";
        return 1;
    }
    std::string file1 = argv[1];
    std::string file2 = argv[2];

    // Read sequences from FASTA files
    std::ifstream fin1(file1);
    std::ifstream fin2(file2);
    if(!fin1.is_open() || !fin2.is_open()) {
        std::cerr << "Error: unable to open input FASTA file(s).\n";
        return 1;
    }
    std::string name1, name2;
    std::string seq1 = "", seq2 = "";
    std::string line;
    // Read first file
    if(std::getline(fin1, line)) {
        if(line.size() > 0 && line[0] == '>') {
            // Extract name (up to first whitespace or end of line after '>')
            size_t pos = 1;
            while(pos < line.size() && !isspace(static_cast<unsigned char>(line[pos]))) {
                name1.push_back(line[pos]);
                pos++;
            }
        }
        // Read sequence lines
        while(std::getline(fin1, line)) {
            if(line.size() > 0 && line[0] == '>') break; // stop if another header (should not happen for single sequence)
            for(char c : line) {
                if(!isspace(static_cast<unsigned char>(c))) {
                    seq1.push_back(c);
                }
            }
        }
    }
    // Read second file
    if(std::getline(fin2, line)) {
        if(line.size() > 0 && line[0] == '>') {
            size_t pos = 1;
            while(pos < line.size() && !isspace(static_cast<unsigned char>(line[pos]))) {
                name2.push_back(line[pos]);
                pos++;
            }
        }
        while(std::getline(fin2, line)) {
            if(line.size() > 0 && line[0] == '>') break;
            for(char c : line) {
                if(!isspace(static_cast<unsigned char>(c))) {
                    seq2.push_back(c);
                }
            }
        }
    }
    fin1.close();
    fin2.close();

    // If names are not provided in FASTA, use filenames instead
    if(name1.empty()) {
        name1 = extractBaseName(file1);
    }
    if(name2.empty()) {
        name2 = extractBaseName(file2);
    }

    // IMPORTANT: Strip any prefix like "BB11001_" from sequence names
    // This change ensures consistent sequence naming in output MSF files
    auto stripPrefix = [](const std::string& name) -> std::string {
        size_t pos = name.find('_');
        if(pos != std::string::npos && pos > 0) {
            return name.substr(pos + 1);
        }
        return name;
    };
    
    // Strip any prefixes from sequence names
    name1 = stripPrefix(name1);
    name2 = stripPrefix(name2);

    // Convert sequences to uppercase (for consistency in matching)
    for(char &c : seq1) {
        c = std::toupper(static_cast<unsigned char>(c));
    }
    for(char &c : seq2) {
        c = std::toupper(static_cast<unsigned char>(c));
    }

    int len1 = seq1.length();
    int len2 = seq2.length();
    if(len1 == 0 || len2 == 0) {
        std::cerr << "Error: one of the sequences is empty.\n";
        return 1;
    }

    // Allocate device memory
    char *d_seq1 = nullptr, *d_seq2 = nullptr;
    int *d_score = nullptr;
    unsigned char *d_dir = nullptr;
    size_t sizeScore = (size_t)(len1+1) * (len2+1) * sizeof(int);
    size_t sizeDir   = (size_t)(len1+1) * (len2+1) * sizeof(unsigned char);
    cudaMalloc((void**)&d_seq1, len1 * sizeof(char));
    cudaMalloc((void**)&d_seq2, len2 * sizeof(char));
    cudaMalloc((void**)&d_score, sizeScore);
    cudaMalloc((void**)&d_dir, sizeDir);
    // Copy sequences to device
    cudaMemcpy(d_seq1, seq1.data(), len1 * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_seq2, seq2.data(), len2 * sizeof(char), cudaMemcpyHostToDevice);
    // Initialize score and direction matrices to 0
    cudaMemset(d_score, 0, sizeScore);
    cudaMemset(d_dir,   0, sizeDir);

    // Scoring scheme (can be adjusted): match = +2, mismatch = -1, gap = -1
    int matchScore = 2;
    int mismatchScore = -1;
    int gapScore = -1;

    // Compute DP matrix anti-diagonal by anti-diagonal
    // Maximum possible diag index = len1 + len2 (when i=len1, j=len2)
    int maxDiag = len1 + len2;
    int threadsPerBlock = 256;
    for(int diag = 2; diag <= maxDiag; ++diag) {
        int start_i = (diag > len2+1) ? (diag - (len2+1) + 1) : 1;
        if(start_i < 1) start_i = 1;
        int end_i = (diag - 1 < len1) ? (diag - 1) : len1;
        if(end_i > len1) end_i = len1;
        if(start_i > len1 || start_i > end_i) continue; // no cells on this diag
        int totalCells = end_i - start_i + 1;
        int blocks = (totalCells + threadsPerBlock - 1) / threadsPerBlock;
        sw_kernel<<<blocks, threadsPerBlock>>>(d_seq1, d_seq2, len1, len2, diag, start_i, end_i, d_score, d_dir, matchScore, mismatchScore, gapScore);
        cudaDeviceSynchronize();
    }

    // Copy score and direction matrices back to host
    std::vector<int> score((len1+1) * (len2+1));
    std::vector<unsigned char> dir((len1+1) * (len2+1));
    cudaMemcpy(score.data(), d_score, sizeScore, cudaMemcpyDeviceToHost);
    cudaMemcpy(dir.data(), d_dir, sizeDir, cudaMemcpyDeviceToHost);

    // Find the cell with maximum score for local alignment endpoint
    int maxScore = 0;
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
    std::string align1 = "";
    std::string align2 = "";
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

    int alignLen = align1.size();
    // Compute checksum for each aligned sequence (GCG checksum)
    auto gcgChecksum = [](const std::string &s) {
        long check = 0;
        for(size_t i = 0; i < s.size(); ++i) {
            // use (i % 57) + 1 multiplier
            int c = std::toupper(static_cast<unsigned char>(s[i]));
            check += ((int)((i % 57) + 1) * c);
        }
        return (int)(check % 10000);
    };
    int check1 = gcgChecksum(align1);
    int check2 = gcgChecksum(align2);
    int globalCheck = ( (long)check1 + (long)check2 ) % 10000;

    // Determine sequence type for MSF (P for protein, N for nucleic acid) by scanning sequence letters
    char typeChar = 'P';
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
    if(maybeDNA) typeChar = 'N';

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

    // Free device memory
    cudaFree(d_seq1);
    cudaFree(d_seq2);
    cudaFree(d_score);
    cudaFree(d_dir);
    return 0;
}