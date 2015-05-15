
# load modules
module load gcc jdk/1.8 tassel/5
# to export hdf5 file to Hapmap format
run_pipeline.pl -Xmx64g -fork1 -h5 ZeaGBSv27_publicSamples_imputedV5_AGPv2-150114.h5 \\
-export -exportType Hapmap -runfork1

# sort the GBS2.7 - this takes a while (1 hour)
run_pipeline.pl -Xmx64g -SortGenotypeFilePlugin -inputFile \\
ZeaGBSv27_publicSamples_imputedV5_AGPv2-150114.hmp.txt -outputFile ZeaGBSv27_sorted \\
-fileType Hapmap

# for GBS2.7, the "keep_list_NAM_children.txt" is just a one column list to keep
run_pipeline.pl -Xmx64g -fork1 -h ZeaGBSv27_sorted.hmp.txt -includeTaxaInfile \\
~/Documents/Github/Misc/data/Taxa_ames282_288.txt -export ZeaGBSv27_Ames282 \\
-exportType Hapmap -runfork1
