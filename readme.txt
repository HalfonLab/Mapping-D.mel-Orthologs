parameters/required files for OM_mappingFlyOrthologsToSCRMshawPredictions
-sp1id speciesX_finalmapped.gff (mapping file between IDs e.g gene-LOC12133 XP214134)
-so scrmshawOutput_speciesX.bed 
-mD DMEL_PROTEIN.fs.maptxt (orthologer output in Rawdata directory)
-mX speciesX_PROTEIN.fs.maptxt  (orthologer output in Rawdata directory)
-og mydata.og_map_speciesX  (orthologer output in Cluster directory)
-ft GCF_000001215.4_Release_6_plus_ISO1_MT_feature_table.txt (file from Github)



./OM_mappingFlyOrthologsToSCRMshawPredictions.py -sp1id species_X_finalMapped.gff -ft GCF_000001215.4_Release_6_plus_ISO1_MT_feature_table.txt -so scrmshawOutput_peaksCalled_allTsets_allMethods_speciesX.bed -mD DMEL_PROTEIN.fs.maptxt -mX speciesX_PROTEIN.fs.maptxt -og mydata.og_map_speciesX


files for CGstoSymbols_final.py
Enter the path to the first file (fbgn_annotation_ID_fb_2023_04.tsv): fbgn_annotation_ID_fb_2023_04.tsv (file from github)
Enter the path to the SCRM file containing CGs:SO_scrmshawOutput_peaksCalled_MedianPointAmpCurve_oldNew48sets_allMethods_A.aeg.bed ( SCRMshaw outout file, output of OM_mappingFlyOrthologsToSCRMshawPredictions.py)
Enter the desired output file name with symbols:SO_scrmshawOutput_peaksCalled_MedianPointAmpCurve_oldNew48sets_allMethods_A.aeg_final.bed (final output name)
