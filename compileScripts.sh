#!/bin/bash
# build_sw.sh - Build script for Smith-Waterman CPU and GPU implementations

# Set colors for terminal output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Custom filenames
cpp_file="cpuSmithWaterman.cpp"
cuda_file="smithWaterman.cu"
cpu_binary="cpuSmithWaterman"
gpu_binary="smithWaterman"

# Compiler options
cpp_compiler="g++"
cuda_compiler="nvcc"
cpp_flags="-O3 -std=c++11"
cuda_flags="-O3 -std=c++11"

# Print header
echo -e "${BLUE}===== Smith-Waterman Implementation Builder =====${NC}"
echo -e "This script compiles both CPU and GPU implementations."
echo ""

# Check if files exist
if [ ! -f "$cpp_file" ]; then
    echo -e "${RED}Error: $cpp_file not found in current directory${NC}"
    exit 1
fi

if [ ! -f "$cuda_file" ]; then
    echo -e "${RED}Error: $cuda_file not found in current directory${NC}"
    exit 1
fi

# Clean previous builds if they exist
if [ -f "$cpu_binary" ]; then
    echo -e "${YELLOW}Removing previous CPU binary...${NC}"
    rm "$cpu_binary"
fi

if [ -f "$gpu_binary" ]; then
    echo -e "${YELLOW}Removing previous GPU binary...${NC}"
    rm "$gpu_binary"
fi

# Check for G++ compiler
if ! command -v $cpp_compiler &> /dev/null; then
    echo -e "${RED}Error: $cpp_compiler compiler not found.${NC}"
    echo -e "Please install g++ to compile the CPU version."
    exit 1
fi

# Check for NVCC compiler
if ! command -v $cuda_compiler &> /dev/null; then
    echo -e "${RED}Error: $cuda_compiler compiler not found.${NC}"
    echo -e "Please ensure CUDA toolkit is installed and in your PATH."
    exit 1
fi

# Build CPU version
echo -e "${BLUE}Building CPU implementation...${NC}"
$cpp_compiler $cpp_flags -o $cpu_binary $cpp_file

if [ $? -eq 0 ]; then
    echo -e "${GREEN}CPU build successful: $cpu_binary${NC}"
else
    echo -e "${RED}CPU build failed${NC}"
    exit 1
fi

# Build GPU version
echo -e "${BLUE}Building GPU implementation...${NC}"
$cuda_compiler $cuda_flags -o $gpu_binary $cuda_file

if [ $? -eq 0 ]; then
    echo -e "${GREEN}GPU build successful: $gpu_binary${NC}"
else
    echo -e "${RED}GPU build failed${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}Both implementations built successfully!${NC}"
echo -e "CPU binary: ${YELLOW}$cpu_binary${NC}"
echo -e "GPU binary: ${YELLOW}$gpu_binary${NC}"
echo ""
echo -e "${BLUE}Usage:${NC}"
echo -e "./$cpu_binary <seq1.fasta> <seq2.fasta>"
echo -e "./$gpu_binary <seq1.fasta> <seq2.fasta>"
echo ""
echo -e "${BLUE}For benchmarking:${NC}"
echo -e "time ./$cpu_binary <seq1.fasta> <seq2.fasta> > cpu_result.txt"
echo -e "time ./$gpu_binary <seq1.fasta> <seq2.fasta> > gpu_result.txt"
echo -e "diff cpu_result.txt gpu_result.txt  # Should show no differences"
echo ""