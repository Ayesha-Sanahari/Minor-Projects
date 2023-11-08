#!/bin/bash

files=$(ls *.fasta)
for i in $files
	do
		pattern=$(grep -B 20 "WGKWVAEIR" $i)
		echo "$pattern" >> pattern_contain.txt		
	done
header=$(grep -e ">" pattern_contain.txt)
echo "$header" >> AP2_basic_headers.txt
