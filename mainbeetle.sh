#!/bin/bash

cd /home/lauralaa/Example/
python /home/lauralaa/Scripts_beetle/get_orthologues.py orthologues.json

mkdir Flies

find Beetles -name 'msa*.fasta' -print0 | xargs -P 5 -0 -I % python /home/lauralaa/Scripts_beetle/get_filtered_fly_data.py orthologues.json % Flies

mkdir Excluded
mv Beetles/* Excluded/

# to make sure sets have matching genes
find Flies -name '*.fasta' -print0 | xargs -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; echo "$base"' -- {} -print0 | sed -e s/[^0-9]//g | xargs -I % find Excluded -name '*[a-z]%.*' | xargs -I % mv % Beetles

# run SLR on the flies, use Ensembl trees
find Flies -name 'fly*.fasta' -print0 | xargs -P 40 -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; python /home/lauralaa/Scripts_beetle/run_slr.py "$1" n > "Flies/${name}.log"' -- {}

# run SLR on the beetles, create trees with RAxML
find Beetles -name 'msa*.fasta' -print0 | xargs -P 40 -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; python /home/lauralaa/Scripts_beetle/run_slr.py "$1" y > "Beetles/${name}.log"' -- {}

# again to make sure sets have matching genes
mkdir Slr_Beetles
mkdir Slr_Flies
find Flies -name '*.slr' -print0 | xargs -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; echo "$base"' -- {} -print0 | sed -e s/[^0-9]//g | xargs -I % find Beetles -name '*[a-z]%.slr*' | xargs -I % mv % Slr_Beetles
find Slr_Beetles -name '*.slr' -print0 | xargs -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; echo "$base"' -- {} -print0 | sed -e s/[^0-9]//g | xargs -I % find Flies -name '*[a-z]%.slr*' | xargs -I % mv % Slr_Flies

# analyze the data
python /home/lauralaa/Scripts_beetle/slr_analysis.py Slr_Beetles beetles
python /home/lauralaa/Scripts_beetle/slr_analysis.py Slr_Flies flies

echo "Graphs in beetles.pdf, flies.pdf"

pwd
