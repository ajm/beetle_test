START

Folder "Example" has the data ready for a test run
Folder "Done_Example" is what you should get after running the scripts on the example data

Folder "Scripts_beetle" contains:
 get_orthologues.py
 get_filtered_flydata.py
 rmgaps.py
 rmtribolium.py
 run_slr_paml.py
 paml_analysis.py
 slr_analysis.py (requires plot.R in the same folder or update the path in the script)
 plot.R

Folder "Beetles" (in the Example folder) has the beetledata:
 msa91.fasta
 msa1955.fasta
 msa2910.fasta
 msa5656.fasta

mainbeetle.sh </path/to/directory/Example/> runs everything, just update paths


RESULTS

 beetles.results
    (SLR omega values of the beetles <= 2)
 flies.results
    (SLR omega values of the flies <=2)
 beetles_means.results
    (SLR omega means of every beetle msa, no filtering)
 flies_means.results
    (SLR omega means of every fly msa, no filtering)
 beetles_fltrd_means.results
    (SLR omega means of every beetle msa, skipping omegas that are over 2)
 flies_fltrd_means.results
    (SLR omega means of every fly msa, skipping omegas that are over 2)

 beetles.pdf
    (3 Slr graphs and 1 Paml graph of beetle results)
 flies.pdf
    (3 Slr graphs and 1 Paml graph of fly results)

Other folders/files after the run is complete:
 Beetles (all the beetle data, info about the SLR run (.log, .stats))
 Flies (fly data, trees from Ensembl, info about the Slr run (.log, .stats))
 Slr_Beetles (the Slr files for the successful beetle Slrs that have a fly counterpart)
 Slr_Flies (the Slr files for the successful fly Slrs that have a beetle counterpart)
 Paml_Beetles (Paml results for beetles)
 Paml_Flies (Paml results for flies)
 Excluded (skipped msa files, because not enough fly or beetle species etc)
 Beetlewithtribolium (original beetle .fasta files before removing tribolium)
 Flyensemblalign (flies aligned by Ensembl)
 Flygapless (fly MSA gaps removed)


