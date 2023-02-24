## script for Hawaiian Lobeliads 

# Trimming 

First, raw data are within the subdirectories of subdirectories, so need to move file up with the following code: 

```
#!/bin/bash
for dir in */; do 
	cd $dir
	for dir2 in */; do 
		cd $dir2 
		mv *.fastq.gz ../
		cd ../../
	done
done
```
Then, some samples contain multiple files and some samples contain regular R1 and R2 files 
Thus, need to "cat" files first. 
A file named "cat_reads" is created first and then the cated files were 

```
#!/bin/bash
# mkdir /Volumes/givnish/Lobeliads_copy_Bing/cat_reads
for dir in */; do 
	cd $dir
	if [ "$(ls -b | wc -l)" -gt 3 ]; then
		cat *_R1_*.fastq.gz > /Volumes/givnish/Lobeliads_copy_Bing/cat_reads/"$(basename "$dir")"_R1.fastq.gz
		cat *_R2_*.fastq.gz > /Volumes/givnish/Lobeliads_copy_Bing/cat_reads/"$(basename "$dir")"_R2.fastq.gz
	else
		mv *_R1_001.fastq.gz /Volumes/givnish/Lobeliads_copy_Bing/cat_reads/"$(basename "$dir")"_R1.fastq.gz
		mv *_R2_001.fastq.gz /Volumes/givnish/Lobeliads_copy_Bing/cat_reads/"$(basename "$dir")"_R2.fastq.gz
	fi
	cd ../
done
```
Then, use FastP to in botany server 
Here, quality score below 20 and read length below 100 were trimmed out, based on the recommendation from Jacob. 
Mostly, more than 90% bases left, so I think the threshold is not too harsh. 

```
# mkdir /staging/bli283/cleaned_reads
for file in /staging/bli283/cat_reads/*R1.fastq.gz
do
	name=`basename $file _R1.fastq.gz`
	echo "Running fastp on $name"
	forward=$name"_R1.fastq.gz"
	reverse=$name"_R2.fastq.gz"
	fastp -i /mnt/researchdrive/givnish/Lobeliads_copy_Bing/cat_reads/$forward -o /mnt/researchdrive/givnish/Lobeliads_copy_Bing/cleaned_reads/$name"_cleaned_R1_.fastq.gz" -I /mnt/researchdrive/givnish/Lobeliads_copy_Bing/cat_reads/$reverse -O /mnt/researchdrive/givnish/Lobeliads_copy_Bing/cleaned_reads/$name"_cleaned_R2_.fastq.gz" -z 4 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 20 --trim_poly_g --length_required 100 --thread 8
done

```

# Annotate the plastome coded region

1. downloaded the reference plastome from genebank (NC_035355)
2. Uploaded the reference plastome to Chloe.plastid.org and annotated it. Then, downloaded .gbff file (Also downloaded .gff3 for reference) 
3. In bbedit, use "replace all" and "extract" function to create a range.txt file (Cl_f_ref_plastome:56142-57614 ...)
4. Editted the reference plastome file to make sure the first line can match with the title of range.txt (Cl_f_ref_plastome)
5. Used the code (SamTools) to only extract coded region from the reference fasta 

```
#!/bin/bash
/home/bli283/SamTools/samtools-1.15.1/samtools faidx Cl_f_ref_plastome.fasta --region-file range.txt > Cl_f_plastome_protein_coding_genes.fasta

```
The output Cl_f_plastome_protein_coding_genes.fasta can then be used as the reference file for HybPiper 

To run hybpiper, need a namelist.txt file 
```
cd /Volumes/givnish/Lobeliads_copy_Bing/cat_reads/
basename *R1.fastq.gz > namelist.txt 
# Then, use bbedit to replace all R1.fastq.gz 

```
Then, I run the following code in the botany server to assemble the exons in plastome. It is noted that for some reasons that I don't understand, fastq.gz extension does not work, so I unzip all fastq.gz files into fastq files before running hybpiper. 

```
per#run HybPaper

# export PATH=/home/bli283/SamTools/samtools-1.15.1/:$PATH
# Always has a SamTools error, so before running the actual script, run the above code fisrt 

while read name; 
do /home/bli283/HybPiper/reads_first.py -b Cl_f_plastome_protein_coding_genes.fasta -r /mnt/researchdrive/givnish/Lobeliads_copy_Bing/cleaned_reads/$name*.fastq --prefix $name --bwa --cpu 6
done < namelist.txt

#retrieve sequences
python /home/bli283/HybPiper/get_seq_lengths.py Cl_f_plastome_protein_coding_genes.fasta namelist.txt dna > test_seq_lengths.txt
# If the script has python in the begginning, the script does not need to be excutable, because it pathyon will make do it and make it as an excutable. 

#stasts
python /home/bli283/HybPiper/hybpiper_stats.py test_seq_lengths.txt namelist.txt > test_stats.txt

#retrieve sequences
python /home/bli283/HybPiper/retrieve_sequences.py Cl_f_plastome_protein_coding_genes.fasta . dna

#cleanup script
python /home/bli283/HybPiper/cleanup.py *

while read name; 
do python /home/bli283/HybPiper/cleanup.py $name
done < namelist.txt
```
Then, I got .FNA files for each loci, and I used the following code to change the extension from .FNA to .fasta 
```
#!/bin/bash
for f in *.FNA; do 
    mv -- "$f" "${f%.FNA}.fasta"
done
```
Then, I did alignment by MAFFT and used trimal to trim the alignment 
```
# doing alignment 
#!/bin/bash
for file in *.fasta
do
	name=`basename $file .fasta`
	mafft --maxiterate 5000 --auto --adjustdirectionaccurately --thread 8 --op 3 --leavegappyregion $file > ~/MAFFT_Bro_baits/$name.MAFFT.fasta
done

# Trim alignment 
#!/bin/bash
for file in *.fasta
do
	name=`basename $file MAFFT.fasta`
	~/miniconda3/bin/trimal -in $file -out trimed_alignment/$name.trimAL.fasta -automated1
done	
```
Then, I used SequenceMatrix to concatenate the files and got a nexus partition file to be used for tree building. 

The partition.nex file is specified when building the IQ-tree.  

I used IQ-tree to build the tree. 

```
iqtree -s Lobeliad_alignment_183.phy -spp partition.nex -m MFP+MERGE -bb 1000 --bsam GENESITE -nt 4  
```
Then, I run the following R-code to change the label tip from sample ID to taxon name 
```
# R script 
library(ape)

setwd("~/Desktop/Lobeliad_changed_tipname")
data <- read.csv("Hawaiian_lobeliads_samples_2021.csv")
tree <- read.tree("partition.nex.treefile")

# Change tree tip 
names(data)
data[[1]] <- sub("\\ ","",data[[1]])
tree$tip.label <- data[[9]][match(tree$tip.label,data[[1]])]
plot(tree)
write.tree(tree, "iqtree_plastome_exons")
```


 



