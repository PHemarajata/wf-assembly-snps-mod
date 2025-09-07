#!/bin/bash

# Bash functions for wf-assembly-snps pipeline

# Function to print messages with timestamp
msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to verify minimum file size
verify_minimum_file_size() {
    local file="$1"
    local description="$2"
    local min_size="$3"
    
    if [[ ! -f "$file" ]]; then
        msg "ERROR: File '$file' does not exist"
        return 1
    fi
    
    # Convert size to bytes
    local min_bytes
    case "$min_size" in
        *k|*K) min_bytes=$((${min_size%[kK]} * 1024)) ;;
        *m|*M) min_bytes=$((${min_size%[mM]} * 1024 * 1024)) ;;
        *g|*G) min_bytes=$((${min_size%[gG]} * 1024 * 1024 * 1024)) ;;
        *) min_bytes="$min_size" ;;
    esac
    
    local actual_size
    actual_size=$(stat -c%s "$file" 2>/dev/null || echo 0)
    
    if [[ "$actual_size" -ge "$min_bytes" ]]; then
        msg "INFO: $description file '$file' passed size check ($actual_size >= $min_bytes bytes)"
        return 0
    else
        msg "WARN: $description file '$file' failed size check ($actual_size < $min_bytes bytes)"
        return 1
    fi
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to get file extension
get_extension() {
    local filename="$1"
    echo "${filename##*.}"
}

# Function to get basename without extension
get_basename() {
    local filename="$1"
    local basename="${filename##*/}"
    echo "${basename%.*}"
}