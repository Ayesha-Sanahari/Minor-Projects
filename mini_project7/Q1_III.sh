#!/bin/bash

files=$(ls *.fasta)
for i in $files
	do
		seq=$(grep -E -r -l "WGKWV|AAEIR" $i)
		echo "$seq" >> AP2_advanced_headers.txt
	done

