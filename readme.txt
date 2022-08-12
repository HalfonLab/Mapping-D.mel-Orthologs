parameters:
-sp1id speciesX_finalmapped.gff 
-so scrmshawOutput_speciesX.bed 
-mD DMEL_PROTEIN.fs.maptxt 
-mX speciesX_PROTEIN.fs.maptxt 
-og mydata.og_map_speciesX 
-ft GCF_000001215.4_Release_6_plus_ISO1_MT_feature_table.txt 



./OM_step1_SCRMsToSCRMswithMappedOrthologs.py -sp1id species_X_finalMapped.gff -ft GCF_000001215.4_Release_6_plus_ISO1_MT_feature_table.txt -so scrmshawOutput_peaksCalled_allTsets_allMethods_speciesX.bed -mD DMEL_PROTEIN.fs.maptxt -mX speciesX_PROTEIN.fs.maptxt -og mydata.og_map_speciesX
