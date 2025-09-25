# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Oligo reconstruction with de bruijn graphs
# Sanchita Bhadra 01/2024
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# project dir:
# rmcox@hfogcomp01:/stor/work/Marcotte/project/rmcox/oligos

# =======================================
# general notes
# =======================================

# 4 baits; 3 around ~40nt, 1 blind sequence
# prey library is ~90nt with N15 variable region flanked by 32nt 5' and 43nt 3' 
# N15 sequence interval: [33,42]

# oligo 1 seq:
# AAGACCTACCTCACATGGCCAACACTCGGACAAAAAAAAAA

# oligo 2 seq:
# GACCACCAGCAGCAGCCAGCCGACGCAGGGACAAAAAAAAAA

# oligo 3 seq:
# AGAACTTACATCAACTAAACAACAAATGAACAAAAAAAAAA

# library sequence
# TATTGCGATAGCTGAGAGAGAAGACGCGAGGGNNNNNNNNNNNNNNNGCGAAAACAAAAAACAAAAATAAGAATCCAAGCAGCAGCAACA

# =======================================
# directory set up
# =======================================

cd /stor/work/Marcotte/project/rmcox/oligos
mkdir {data,cmds,logs,results,figures,scripts}
mkdir -p data/raw/{oligo_1,oligo_2,oligo_3,oligo_b}
mkdir -p data/qc/{oligo_1,oligo_2,oligo_3,oligo_b,fastp}
mkdir -p data/merged/{oligo_1,oligo_2,oligo_3,oligo_b}
mkdir -p data/processed/{oligo_1,oligo_2,oligo_3,oligo_b}
mkdir -p results/{oligo_1,oligo_2,oligo_3,oligo_b}/{graphs,counts,enrichment}

# view directory structure
tree data/

# =======================================
# symbolic link raw data
# =======================================

# oligo
ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Oligo1* data/raw/oligo_1
ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Oligo2* data/raw/oligo_2
ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Oligo3* data/raw/oligo_3
ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*OligoBlind* data/raw/oligo_b

# =======================================
# data QC
# [01/2024 - 02/2024]
# =======================================

# >>>>>>>>>>>>>>>>>>>>>>>
# running fastp [01/2024]
# >>>>>>>>>>>>>>>>>>>>>>>

# template command for fastp:
# fastp -i R1.fastq.qz -I R2.fastq.gz -o R1.trimmed.fastq.gz -o R2.trimmed.fastq.gz -p --failed_out data/qc/fastp_failed/

# print file paths to R1 (left) reads
for exp in 1 2 3 b;
	do ls -C data/raw/oligo_${exp}/*R1*
done

# for each oligo exp, write file paths to R1 reads into a new file
mkdir -p data/meta/raw_read_paths/
for x in 1 2 3 b; do ls -C data/raw/oligo_${x}/*_R1_001.fastq.gz > data/meta/raw_read_paths/oligo_${x}_R1.txt; done

# for each oligo exp, write file paths to R2 reads into a new file
for x in 1 2 3 b; do ls -C data/raw/oligo_${x}/*_R2_001.fastq.gz > data/meta/raw_read_paths/oligo_${x}_R2.txt; done

# check to make sure same number of R1 and R2 reads
# for each exp 
wc -l data/meta/raw_read_paths/*txt

# join file paths together into single file
for x in 1 2 3 b; do paste data/meta/raw_read_paths/oligo_${x}_R1.txt data/meta/raw_read_paths/oligo_${x}_R2.txt > data/meta/raw_read_paths/oligo_${x}_R1R2.txt; done

# creating the fastp command loop
mkdir -p data/qc/fastp/

# step 1: loop over files to create names for input and output
for x in 1 2 3 b; do
	while IFS=$'\t' read -r -a y; do
	echo "R1: ${y[0]} | R2: ${y[1]}"
    echo "O1: ${y[0]##*/} | O2: ${y[1]##*/}"
    done < data/meta/raw_read_paths/oligo_${x}_R1R2.txt
done

# step 2: loop over files to create names for input and output AND create names for the merged output
for x in 1 2 3 b; do
	while IFS=$'\t' read -r -a y; do
    	echo "R1: ${y[0]} | R2: ${y[1]}"
    	echo "O1: ${y[0]##*/} | O2: ${y[1]##*/}"
    	b1="$(basename ${y[0]} .fastq.gz)"
    	b2="$(basename ${y[1]} .fastq.gz)"
    	echo "B1: data/qc/oligo_${x}/${b1}.fpqc.fastq.gz | B2: data/qc/oligo_${x}/${b2}.fpqc.fastq.gz"
    	merg="${b1%_L00*}.merged.fpqc.fastq.gz"
    	echo "M: ${merg}"
	done < data/meta/raw_read_paths/oligo_${x}_R1R2.txt
done

# step 3: fill in filenames for the fastp command & then write each fastp command to file
for x in 1 2 3 b; do
	i=0
	while IFS=$'\t' read -r -a y; do
    	o1="$(basename ${y[0]} .fastq.gz)"
    	o2="$(basename ${y[1]} .fastq.gz)"
    	merg="${o1%_L00*}.merged.fpqc.fastq.gz"
    	echo "fastp -i ${y[0]} -I ${y[1]} -o data/qc/oligo_${x}/${o1}.fpqc.fastq.gz -O data/qc/oligo_${x}/${o2}.fpqc.fastq.gz --failed_out data/qc/fastp/failed/${o1%_L00*}.failed.txt --merge --merged_out data/merged/oligo_${x}/${merg} --length_limit 100 --json data/qc/fastp/oligo_${x}_${o1%_L00*}.json --html data/qc/fastp/oligo_${x}_${o1%_L00*}.html"
	done < data/meta/raw_read_paths/oligo_${x}_R1R2.txt
done > cmds/fastp_oligos.cmds

# run 12 fastp commands at a time
cat cmds/fastp_oligos.cmds | parallel -j12

# >>>>>>>>>>>>>>>>>>>>>>>
# running multiqc [02/2024]
# >>>>>>>>>>>>>>>>>>>>>>>

# use multiqc to summarize the results of fastp trim & merge

# generate multiqc reports for oligos
for x in 1 2 3 b; do multiqc data/qc/fastp/oligo_${x}* --filename data/qc/multiqc/oligo_${x}_merged --interactive --force; done

# =======================================
# renaming (organizing) merged reads +
# decompressing and converting to FASTA
# [02/2024]
# =======================================

mkdir -p data/meta/merged_read_paths/

# get paths & number replicates for each experiment (oligos)
for x in 1 2 3 b; do
	grep -n '.*' <<< `ls data/merged/oligo_${x}/*Probe*` > data/meta/merged_read_paths/oligo_${x}_bg_file_paths.txt
	grep -n '.*' <<< `ls data/merged/oligo_${x}/*0w*` > data/meta/merged_read_paths/oligo_${x}_0w_file_paths.txt
	grep -n '.*' <<< `ls data/merged/oligo_${x}/*1w*` > data/meta/merged_read_paths/oligo_${x}_1w_file_paths.txt
	grep -n '.*' <<< `ls data/merged/oligo_${x}/*3w*` > data/meta/merged_read_paths/oligo_${x}_3w_file_paths.txt
done

# >>>>>>>>>>>>>>>>>>>>>>>
# background data
# >>>>>>>>>>>>>>>>>>>>>>>
# oligos
for x in 1 2 3 b; do
	while IFS=$':' read -r -a y; do
	rep_num=${y[0]}
	file=${y[1]}
	new_name="oligo_${x}_bg_rep${rep_num}.merged.fasta"
	zcat $file | sed -n '1~4s/^@/>/p;2~4p' > data/processed/oligo_${x}/${new_name}
done < data/meta/merged_read_paths/oligo_${x}_bg_file_paths.txt
done

# >>>>>>>>>>>>>>>>>>>>>>>
# 0w wash data
# >>>>>>>>>>>>>>>>>>>>>>>
# oligos
for x in 1 2 3 b; do
	while IFS=$':' read -r -a y; do
	rep_num=${y[0]}
	file=${y[1]}
	new_name="oligo_${x}_0w_rep${rep_num}.merged.fasta"
	zcat $file | sed -n '1~4s/^@/>/p;2~4p' > data/processed/oligo_${x}/${new_name}
done < data/meta/merged_read_paths/oligo_${x}_0w_file_paths.txt
done

# >>>>>>>>>>>>>>>>>>>>>>>
# 1w wash data
# >>>>>>>>>>>>>>>>>>>>>>>
# oligos
for x in 1 2 3 b; do
	while IFS=$':' read -r -a y; do
	rep_num=${y[0]}
	file=${y[1]}
	new_name="oligo_${x}_1w_rep${rep_num}.merged.fasta"
	zcat $file | sed -n '1~4s/^@/>/p;2~4p' > data/processed/oligo_${x}/${new_name}
done < data/meta/merged_read_paths/oligo_${x}_1w_file_paths.txt
done

# >>>>>>>>>>>>>>>>>>>>>>>
# 3w wash data
# >>>>>>>>>>>>>>>>>>>>>>>
# oligos
for x in 1 2 3 b; do
	while IFS=$':' read -r -a y; do
	rep_num=${y[0]}
	file=${y[1]}
	new_name="oligo_${x}_3w_rep${rep_num}.merged.fasta"
	zcat $file | sed -n '1~4s/^@/>/p;2~4p' > data/processed/oligo_${x}/${new_name}
done < data/meta/merged_read_paths/oligo_${x}_3w_file_paths.txt
done

# =======================================
# generate de Bruijn graphs for oligos
# [02/2024 - 03/2024]
# =======================================

# always assign this "script_dir" var before running repo code
script_dir="/stor/work/Marcotte/project/rmcox/deBruijn_graphs/scripts"
	
# >>>>>>>>>>>>>>>>>>>>>>>
# generate all kmers
# >>>>>>>>>>>>>>>>>>>>>>>

# low memory, quick (minutes to hours)
for k in {5..15}; do
for dir in oligo_1 oligo_2 oligo_3 oligo_b; do
	for infile in `ls data/processed/${dir}`; do
	outfile=${infile%.merged*}_${k}mers.csv
	echo "python3 ${script_dir}/generate_deBruijn_kmers.py \
	--input_fasta data/processed/${dir}/${infile} \
	--outfile results/${dir}/graphs/${outfile} \
	--kmer_size ${k}"
	done
	done
done > cmds/generate_dbgs.sh
cat cmds/generate_dbgs.sh | parallel -j36

# >>>>>>>>>>>>>>>>>>>>>>>
# count kmers
# >>>>>>>>>>>>>>>>>>>>>>>

# memory intensive, takes awhile (hours to days)
for k in {5..15}; do
for dir in oligo_1 oligo_2 oligo_3 oligo_b; do
	for infile in `ls results/${dir}/graphs/*${k}mers.csv`; do
	name="$(basename $infile .csv)"
	outfile="${name}_counts.csv"
	echo "python3 ${script_dir}/collapse_kmer_uniques.py \
	--input_file ${infile} \
	--outfile results/${dir}/counts/${outfile} \
	| tee -a logs/${dir}_${name}_counts.log"
	done
	done
done > cmds/generate_kmer_counters.sh
cat cmds/generate_kmer_counters.sh | parallel -j3

# =======================================
# compute enrichment sans library sequences
# [03/2024]
# =======================================

# make dir for filtered k-mer set
# make dirs for combined results
for dir in oligo_1 oligo_2 oligo_3 oligo_b; do
	mkdir results/${dir}/counts_combined_nolib/
done

# # >>>>>>>>>>>>>>>>>>>>>>>
# label count data & remove lib constant regions
# >>>>>>>>>>>>>>>>>>>>>>>
# template/test
script_dir="/stor/work/Marcotte/project/rmcox/deBruijn_graphs/scripts"
python3 ${script_dir}/label_kmers.py \
--results_file results/oligo_1/counts_combined/oligo_1_10mers_combined_10c_wide.csv --lib_file data/bait_lib_kmers/library_10mers.csv \
--bait_file data/bait_lib_kmers/oligo_1_comp_10mers.csv \
--outfile results/oligo_1/counts_combined_nolib/oligo_1_10mers_combined_10c_wide_nolib.csv \
--remove_lib --drop_label_col

# command generator lib
for exp in oligo_1 oligo_2 oligo_3 oligo_b; do
	for k in {5..15}; do
	for file in `ls results/${exp}/counts_combined/*_${k}mers_*wide*`; do
	fname="$(basename ${file} .csv)"
	echo "python3 ${script_dir}/label_kmers.py \
	--results_file $file \
	--lib_file data/bait_lib_kmers/library_${k}mers.csv \
	--bait_file data/bait_lib_kmers/${exp}_comp_${k}mers.csv \
	--outfile results/${exp}/counts_combined_nolib/${fname}_nolib.csv \
	--remove_lib --drop_label_col"
done
done
done > cmds/remove_lib_kmers.cmds
cat cmds/remove_lib_kmers.cmds | parallel -j16

# >>>>>>>>>>>>>>>>>>>>>>>
# re-run degust on data (sans lib constant regions) to get enriched k-mers
# >>>>>>>>>>>>>>>>>>>>>>>

# make dirs for new data
for dir in oligo_1 oligo_2 oligo_3 oligo_b; do
	mkdir results/${dir}/enrichment_nolib/
done

# command loop (run test for all washes combined)
for k in {5..15}; do
	for exp in oligo_1 oligo_2 oligo_3 oligo_b; do
	echo "Rscript ${script_dir}/run_degust.R \
	--counts_file results/${exp}/counts_combined_nolib/${exp}_${k}mers_combined_10c_wide_nolib.csv \
	--exp 0w,1w,3w \
	--output_file results/${exp}/enrichment_nolib/${exp}_${k}mers_10c_all_washes_enriched_nolib.csv | tee -a logs/${exp}_${k}mers_all_wash_enrichment_nolib.log"
done
done > cmds/calc_enrich_all_wash_oligos_nolib.cmds
cat cmds/calc_enrich_all_wash_oligos_nolib.cmds | parallel -j16

# >>>>>>>>>>>>>>>>>>>>>>>
# label data for plots
# >>>>>>>>>>>>>>>>>>>>>>>

# command generator loop
for exp in oligo_1 oligo_2 oligo_3 oligo_b; do
	for k in {5..15}; do
	for file in `ls results/${exp}/enrichment_nolib/*_${k}mers_*`; do
	fname="$(basename ${file} .csv)"
	echo "python3 ${script_dir}/label_kmers.py --results_file $file --lib_file data/bait_lib_kmers/library_${k}mers.csv --bait_file data/bait_lib_kmers/${exp}_comp_${k}mers.csv --outfile results/${exp}/labeled_nolib/${fname}_labeled.csv"
done
done
done > cmds/label_enriched_kmers_nolib.cmds
cat cmds/label_enriched_kmers_nolib.cmds | parallel -j32

# analysis & plot code in:
# scripts/analyze_data.R

# =======================================
# generate consensus sequences
# [04/2024]
# =======================================

script_dir="/stor/work/Marcotte/project/rmcox/deBruijn_graphs/scripts"

# make dirs for new data
for dir in oligo_1 oligo_2 oligo_3 oligo_b; do
	mkdir -p results/${dir}/consensus_seqs_nolib/
done

# consensus seq command loop
for exp in oligo_1 oligo_2 oligo_3 oligo_b; do
	for k in {5..15}; do
	echo "python3 ${script_dir}/calc_consensus_seq.py --input_file results/${exp}/consensus_seqs_nolib/${exp}_${k}mers_weights.csv --calc_n 100 --edges"
done
done > cmds/top_kmer_weights_nolib.cmds
cat cmds/top_kmer_weights_nolib.cmds | parallel -j8

# =======================================
# label normalized counts
# [06/2025]
# =======================================

# make combined lib/bait files
awk 'NR==1 || FNR>1' data/bait_lib_kmers/oligo_1_comp*.csv > data/bait_lib_kmers/oligo_1_comp_all_kmers.csv
awk 'NR==1 || FNR>1' data/bait_lib_kmers/oligo_2_comp*.csv > data/bait_lib_kmers/oligo_2_comp_all_kmers.csv
awk 'NR==1 || FNR>1' data/bait_lib_kmers/oligo_3_comp*.csv > data/bait_lib_kmers/oligo_3_comp_all_kmers.csv
awk 'NR==1 || FNR>1' data/bait_lib_kmers/oligo_b_comp*.csv > data/bait_lib_kmers/oligo_b_comp_all_kmers.csv
awk 'NR==1 || FNR>1' data/bait_lib_kmers/library*.csv > data/bait_lib_kmers/library_all_kmers.csv


script_dir="/stor/work/Marcotte/project/rmcox/deBruijn_graphs/scripts"
python3 ${script_dir}/label_kmers.py \
--results_file results/normalized_counts_20c/oligo_1_20c_cpms_all_kmers.csv \
 --lib_file data/bait_lib_kmers/library_all_kmers.csv \
 --bait_file data/bait_lib_kmers/oligo_1_comp_all_kmers.csv \
 --outfile results/normalized_counts_20c/oligo_1_20c_cpms_all_kmers_labeled.csv \
 --remove_lib --drop_label_col

for x in 1 2 3 b; do
	echo "python3 ${script_dir}/label_kmers.py \
	--results_file results/normalized_counts_20c/oligo_${x}_20c_cpms_all_kmers.csv \
	--lib_file data/bait_lib_kmers/library_all_kmers.csv \
	--bait_file data/bait_lib_kmers/oligo_${x}_comp_all_kmers.csv \
	--outfile results/normalized_counts_20c/oligo_${x}_20c_cpms_all_kmers_labeled.csv \
	--remove_lib"
done > cmds/label_normalized_counts.cmds
cat cmds/label_normalized_counts.cmds | parallel -j4


# ==================================================================
# -------------------------------end--------------------------------
# ==================================================================
















