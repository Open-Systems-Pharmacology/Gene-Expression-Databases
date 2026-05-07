#!/bin/bash

## Change to the PK-Sim DBs directory
cd "PK-Sim DBs/"
# Find all .expressionDB files recursively and compress each to .7z
#find . -type f -name "*ADME_ONLY*.expressionDB" | while read dbfile; do
#  # Remove leading ./ from filename for output
#  outname="${dbfile#./}.7z"
#  # Compress using 7z (overwrite if exists)
#  7z a -t7z "$outname" "$dbfile"
#done

find . -type f -name "*ADME_ONLY*.expressionDB" | while read dbfile; do
  outname="${dbfile#./}.tar.gz"
  tar -czf "$outname" "$dbfile"
done