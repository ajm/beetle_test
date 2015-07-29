START

Folder "Example" has the data ready for a test run
Folder "Done_Example" is what you should get after running the scripts on the example data

Folder "Scripts_beetle" contains:
 get_orthologues.py
 get_filtered_flydata.py
 run_slr.py
 slr_analysis.py (requires plot.R in the same folder or update the path in the script)
 plot.R

Folder "Beetles" (in the Example folder) has the beetledata:
 msa97.fasta
 msa4291.fasta
 msa431.fasta
 msa1710.fasta
 msa54.fasta
 msa2055.fasta

mainbeetle.sh runs everything, just update paths



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
    (graphs of the 3 beetle result files mentioned above, all values don't show)
 flies.pdf
    (graphs of the 3 fly results files mentioned above)

Other folders/files after the run is complete:
 Beetles (all the beetle data, generated trees etc, info about the SLR run (.log, .stats))
 Flies (fly data, trees from Ensembl, info about the SLR run (.log, .stats))
 Slr_Beetles (the SLR files for the successful beetle SLRs that have a fly counterpart)
 Slr_Flies (the SLR files for the successful fly SLRs that have a beetle counterpart)
 Excluded (skipped msa files, because not enough fly or beetle species etc)



In case of any questions or alarmingly bad code, please contact me.
