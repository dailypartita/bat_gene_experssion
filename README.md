# Bat gene expression matrix

> This instruction only for the Egyptian fruit bat (Rousettus aegyptiacus)

1. Build hisat2 index: `1-work-build_all_index.sh`
2. QC, mapping and then calculate FPKM: `2-work-mapping_each_sample.sh`
3. Count the FPKM of all samples: `3-work-make_summary.sh`
4. Usage
   1. Generate `2-work-mapping_each_sample.sh` with the sample_name, reads1 path, reads2 path (FASTQ) of all samples and use the parameters given by `ignition.sh` (8 cores, 2G ARM), the script automatically clears all the intermediate large files generated in temp.
   2. Run `3-work-make_summary.sh` locally to count all the results, and save as a CSV file: `output/merged_htc.csv`.
5. Running Time: Each sample takes less than 20 minutes.
