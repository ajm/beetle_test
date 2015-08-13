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
 msa25.fasta
 msa2910.fasta
 msa42.fasta
 msa511.fasta
 msa561.fasta
 msa91.fasta
 msa929.fasta

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
 Slr_Beetles (accepted beetle .slr files)
 Slr_Flies (accepted fly .slr files)
 Paml_Beetles (accepted beetle .paml files)
 Paml_Flies (accepted fly .paml files)
 Excluded (skipped msa files, because not enough fly or beetle species etc)
 Beetlewithtribolium (original beetle .fasta files before removing tribolium)
 Flyensemblalign (flies aligned by Ensembl)
 Flygapless (fly MSA gaps removed)


