#!/bin/bash

num=$(cat *.fasta | grep -o '>' | wc -l )
echo "total number of files: $num"

