#!/bin/bash

cd $1
python /home/lauralaa/Scripts_beetle/get_orthologues.py orthologues.json

mkdir Flies

find Beetles -name 'msa*.fasta' -print0 | xargs -P 5 -0 -I % python /home/lauralaa/Scripts_beetle/get_filtered_fly_data.py orthologues.json % Flies

mkdir Excluded
mv Beetles/* Excluded/

# to make sure sets have matching genes
find Flies -name '*.fasta' -print0 | xargs -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; echo "$base"' -- {} -print0 | sed -e s/[^0-9]//g | xargs -I % find Excluded -name '*[a-z]%.*' | xargs -I % mv % Beetles

# remove tribolium from the beetle files
mkdir Beetlewithtribolium
mv Beetles/*.fasta Beetlewithtribolium
find Beetlewithtribolium -name 'msa*.fasta' -print0 | xargs -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; python /home/lauralaa/Scripts_beetle/rmtribolium.py "$1" "Beetles/${name}.fasta"' -- {}

# run SLR on the flies, use Ensembl trees
find Flies -name 'fly*.fasta' -print0 | xargs -P 40 -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; python /home/lauralaa/Scripts_beetle/run_slr_paml.py "$1" n > "Flies/${name}.log"' -- {}

# run SLR on the beetles, create trees with RAxML
find Beetles -name 'msa*.fasta' -print0 | xargs -P 40 -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; python /home/lauralaa/Scripts_beetle/run_slr_paml.py "$1" y > "Beetles/${name}.log"' -- {}

# again to make sure sets have matching genes, create own folders for Slr results
mkdir Slr_Beetles
mkdir Slr_Flies
find Flies -name '*.slr' -print0 | xargs -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; echo "$base"' -- {} -print0 | sed -e s/[^0-9]//g | xargs -I % find Beetles -name '*[a-z]%.slr*' | xargs -I % mv % Slr_Beetles
find Slr_Beetles -name '*.slr' -print0 | xargs -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; echo "$base"' -- {} -print0 | sed -e s/[^0-9]//g | xargs -I % find Flies -name '*[a-z]%.slr*' | xargs -I % mv % Slr_Flies

# to make sure Paml matches with Slr, make own folders for Paml results as well
mkdir Paml_Beetles
mkdir Paml_Flies
find Slr_Beetles -name '*.slr' -print0 | xargs -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; echo "$base"' -- {} -print0 | sed -e s/[^0-9]//g | xargs -I % find Beetles -name '*[a-z]%.paml*' | xargs -I % mv % Paml_Beetles
find Slr_Flies -name '*.slr' -print0 | xargs -0 -I {} sh -c 'base=$(basename $1) ; name=${base%.*} ; ext=${base##*.} ; echo "$base"' -- {} -print0 | sed -e s/[^0-9]//g | xargs -I % find Flies -name '*[a-z]%.paml*' | xargs -I % mv % Paml_Flies

# analyze the data
python /home/lauralaa/Scripts_beetle/paml_analysis.py Paml_Beetles/ beetles_paml.results
python /home/lauralaa/Scripts_beetle/paml_analysis.py Paml_Flies/ flies_paml.results
python /home/lauralaa/Scripts_beetle/slr_analysis.py Slr_Beetles beetles
python /home/lauralaa/Scripts_beetle/slr_analysis.py Slr_Flies flies

echo "Graphs in beetles.pdf, flies.pdf"

pwd
