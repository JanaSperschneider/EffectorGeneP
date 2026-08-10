# `EffectorGeneP`: gene annotation in pathogen genomes

#### What is `EffectorGeneP` for?
Gene annotation is crucial for accurate inference of biological knowledge from genomes. Automated gene annotation methods rely on decades-old methods biased towards model species and conserved genes, whilst manual curation is inefficient and can be error-prone. Incorporation of transcription evidence such as RNA-seq data has vastly improved accurate gene annotation. However, non-canonical genes such as those lacking homologs, those residing in rapidly evolving genomic regions or single-exon genes are still routinely dismissed in annotation pipelines as transcriptional noise. In pathogen genomes, this disproportionately affects the accurate annotation of genes encoding disease-causing effector proteins. 

#### What does `EffectorGeneP` do?
`EffectorGeneP` is a machine learning tool that self-trains on infection transcript data to distinguish coding sequences of non-secreted proteins, secreted proteins and effector proteins from non-coding regions of the genome. `EffectorGeneP` predicts the most likely coding sequence from transcripts and effectively addresses transcript fusions which frequently occur in compact genomes, whilst separating genes from transcriptional noise. In benchmarking, `EffectorGeneP` annotates over 95% of known effectors correctly, while other methods only annotate 15%-78%. `EffectorGeneP` is designed to work for **eukaryotic genomes, particularly fungal genomes**.

#### What do I need for running `EffectorGeneP`?

✅ A fungal (or at least eukaryotic) pathogen genome

✅ Infection transcriptome data (RNA-Seq or Iso-Seq)

<img width="836" height="218" alt="GitHub_WorkFlow_Figure" src="https://github.com/user-attachments/assets/bc66bfe0-477c-4f3d-9db6-0823e34a10fa" />

# Table of contents
* [Installing `EffectorGeneP`](#installing-effectorgenep)
* [Testing `EffectorGeneP` with the toy example](#testing-effectorgenep-with-the-toy-example)
* [Prepare your own input data](#prepare-your-own-input-data)
* [Output files](#output-files)
* [`EffectorGeneP` parameter reference](#effectorgenep-parameter-reference)
* [FAQs](#faqs)
* [Licence](#licence)
* [Contact](#contact)
* [Citation](#citation)

## Installing `EffectorGeneP`

`EffectorGeneP` has been written in Python3 and has the following base requirements:

* Python3, BioPython and SciPy (for general processing)
* Java for running the machine learning classification (WEKA)

To get `EffectorGeneP` to work on your local machine, follow these steps:

**Step 1**. Download the latest release from this github repo (i.e. download EffectorGeneP_1.0.0.zip, *not* the Source code zip/tar.gz files) or alternatively for the latest version you can simply clone the github repo (git clone https://github.com/JanaSperschneider/EffectorGeneP.git).

Unpack `EffectorGeneP` in your desired location if you downloaded a release:
```
unzip EffectorGeneP_1.0.0
```
Then cd into the directory:
```
cd EffectorGeneP_1.0.0
```

**Step 2**. Download the `EffectorGeneP` model files from [here](https://effectorp.csiro.au/effectorgenep.html). **It is best to save and unzip them in the `EffectorGeneP` directory folder ./EffectorGeneP_Models/.** 

:exclamation:*This step is optional but desirable to enable sensitive search for genes encoding secreted proteins*. 

**Optional Step 3**. Obtain and install SignalP 4.1 from [here](https://services.healthtech.dtu.dk/services/SignalP-4.1/) and TMHMM 2.0 from [here](https://services.healthtech.dtu.dk/services/TMHMM-2.0/). **Please do not install newer/older versions of SignalP or TMHMM, it is important to use these for accurate and sensitive effector gene annotation and compatibility of output formats**. 

## Testing `EffectorGeneP` with the toy example
`EffectorGeneP` requires these input files:

* The pathogen genome FASTA file
* A GFF3 file of assembled transcripts in gffread format (see detailed explanation below)
* The `EffectorGeneP` machine learning models for the species, or a related species

To test your installation, run `EffectorGeneP` on the toy example as follows.

First, the input file is the transcript for the MC69 effector gene, have a look at the required gff3 format:
```
head toy_example/MC69_Input.gff3
NC_017854.1     StringTie       transcript      1289343 1290545 1000    -       .       ID=MSTRG.19800.1;geneID=MSTRG.19800
NC_017854.1     StringTie       exon    1289343 1289803 1000    -       .       Parent=MSTRG.19800.1
NC_017854.1     StringTie       exon    1289890 1290545 1000    -       .       Parent=MSTRG.19800.1
```

Then run `EffectorGeneP` as follows:
```
python EffectorGeneP.py -t ./toy_example/MC69_Input.gff3 \
                        -g ./toy_example/NC_017854.1.fasta \
                        -m ./EffectorGeneP_Models/Magnaporthe_oryzae/ \
                        -o out.gff3
```
`EffectorGeneP` will return the GFF3 of annotated genes in the transcript. Have a look at the `EffectorGeneP` GFF3 output file:

```
head out.gff3
##gff-version 3
NC_017854.1     EffectorGeneP   mRNA    1289779 1290029 0.83    -       .       ID=MSTRG.19800.1_ORF1_-;geneID=MSTRG.19800.1
NC_017854.1     EffectorGeneP   CDS     1289779 1289803 .       -       1       ID=cds.MSTRG.19800.1_ORF1_-;Parent=MSTRG.19800.1_ORF1_-
NC_017854.1     EffectorGeneP   CDS     1289890 1290029 .       -       0       ID=cds.MSTRG.19800.1_ORF1_-;Parent=MSTRG.19800.1_ORF1_-
```
In this case, `EffectorGeneP` returns one gene for this transcript. Column 6 in the gff3 file is the `EffectorGeneP` score for encoding an effector gene [0-1], in this case the score is 0.83. Note that if you intend to use the EffectorGeneP scores for ranking, this should only be done using the ML model of the species. 

If you have installed SignalP 4.1 and TMHMM 2.0, you can run `EffectorGeneP` as follows:

```
python EffectorGeneP.py -t ./toy_example/MC69_Input.gff3
                        -g ./toy_example/NC_017854.1.fasta
                        -m ./EffectorGeneP_Models/Magnaporthe_oryzae/
                        -o out.gff3
                        --SIGNALP4 /path/to/signalp4.1/signalp
                        --TMHMM /path/to/tmhmm2.0/bin/tmhmm
```
In addition to the GFF3 file above, `EffectorGeneP` will also return a GFF3 of **annotated genes encoding secreted proteins** in the transcript as predicted by SignalP 4.1 and TMHMM 2.0. Here, the gene found for the transcript encodes a secreted protein, the MC69 effector gene.

```
head out.secretome.gff3
##gff-version 3
NC_017854.1     EffectorGeneP   mRNA    1289779 1290029 0.83    -       .       ID=MSTRG.19800.1_ORF1_-;geneId=MSTRG.19800.1
NC_017854.1     EffectorGeneP   CDS     1289779 1289803 .       -       1       ID=cds.MSTRG.19800.1_ORF1_-;Parent=MSTRG.19800.1_ORF1_-
NC_017854.1     EffectorGeneP   CDS     1289890 1290029 .       -       0       ID=cds.MSTRG.19800.1_ORF1_-;Parent=MSTRG.19800.1_ORF1_-
```

## Prepare your own input data

<img width="3704" height="1148" alt="GitHub_WorkFlow_Figure_Alignment" src="https://github.com/user-attachments/assets/48fdd0b1-c3ee-48d6-9ab4-29a82a6243b5" />

The first step is to gather your RNA-seq infection data and align each replicate/sample to your genome. For example, for gene-dense fungal genomes one could use STAR (https://github.com/alexdobin/STAR) in 2-pass mode as follows:

```
# Align replicate1
STAR --runMode alignReads --genomeDir genome_index/ --outFileNamePrefix Rep1 --readFilesIn ${Rep1_read1} ${Rep2_read2} --alignIntronMin 5 --alignIntronMax 3000 --alignMatesGapMax 3000 --outFilterMultimapNmax 100 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif
samtools index Rep1_Aligned.sortedByCoord.out.bam
# Remove junctions supported by very few reads (e.g. = 2 reads)
cat Rep1_SJ.out.tab | awk '($7 > 2)' | cut -f1-6 | sort | uniq > Rep1_SJ.out.filtered.tab
# Align again with junction information
STAR --runMode alignReads --genomeDir genome_index/ --outFileNamePrefix Rep1_2pass --readFilesIn ${Rep1_read1} ${Rep1_read2} --alignIntronMin 5 --alignIntronMax 3000 --alignMatesGapMax 3000 --outFilterMultimapNmax 100 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --sjdbFileChrStartEnd Rep1_SJ.out.filtered.tab
samtools index Rep1_2pass_Aligned.sortedByCoord.out.bam

# Align replicate2
STAR --runMode alignReads --genomeDir genome_index/ --outFileNamePrefix Rep2 --readFilesIn ${Rep2_read1} ${Rep1_read2} --alignIntronMin 5 --alignIntronMax 3000 --alignMatesGapMax 3000 --outFilterMultimapNmax 100 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif
samtools index Rep2_Aligned.sortedByCoord.out.bam
# Remove junctions supported by very few reads (e.g. = 2 reads)
cat Rep2_SJ.out.tab | awk '($7 > 2)' | cut -f1-6 | sort | uniq > Rep2_SJ.out.filtered.tab
# Align again with junction information
STAR --runMode alignReads --genomeDir genome_index/ --outFileNamePrefix Rep2_2pass --readFilesIn ${Rep2_read1} ${Rep2_read2} --alignIntronMin 5 --alignIntronMax 3000 --alignMatesGapMax 3000 --outFilterMultimapNmax 100 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --sjdbFileChrStartEnd Rep2_SJ.out.filtered.tab
samtools index Rep2_2pass_Aligned.sortedByCoord.out.bam
```

The next step is to assemble transcripts for each replicate **using appropriate strand settings for your data (fr/rf/none)** and then to merge them into a consensus set. Here we use with StringTie (https://github.com/gpertea/stringtie). 
```
# Assemble transcripts from aligned reads in BAM format, this example is --rf stranded
stringtie -s 1 -m 200 Rep1_2pass_Aligned.sortedByCoord.out.bam.bam -o StringTie.Rep1.gtf --rf
stringtie -s 1 -m 200 Rep2_2pass_Aligned.sortedByCoord.out.bam.bam -o StringTie.Rep2.gtf --rf

# Merge into consensus transcript set
stringtie --merge -l MSTRG -o StringTie.gtf StringTie.Rep1.gtf StringTie.Rep2.gtf

# It is essential to convert gtf to gff3 with gffread (https://github.com/gpertea/gffread) for input into EffectorGeneP
gffread StringTie.gtf > EffectorGeneP_input_transcripts.gff3
```

:exclamation:Note: you have to use gffread (https://github.com/gpertea/gffread) to format your gtf/gff3 transcript input file before you input it into `EffectorGeneP`. **If you skip this step, `EffectorGeneP` might not parse your input file correctly.**:exclamation:

Lastly, you can run `EffectorGeneP` on the transcripts as follows:

```
python EffectorGeneP.py [-options] -g <genome fasta file>
                                   -t <transcripts in GFF3 format (gffread formatted)>
                                   -m <path to EffectorGeneP model files>
                                   -o <EffectorGeneP GFF3 output file>
```
:exclamation:Note that only the *Magnaporthe oryzae* model files are included in the Github repository, more models are available [here](https://effectorp.csiro.au/effectorgenep.html).

## Output files

`EffectorGeneP` will use the transcript identifiers to name the genes it finds. `EffectorGeneP` appends a nomenclature to the transcript identifier that indicates if it annotated the longest ORF in the transcript (ORF0) or the second-longest ORF in the transcript (ORF1) etc. as well as the strand of the gene (e.g. ORF1_- or ORF1_+). In column 6 of the gff3 file is the `EffectorGeneP` score for encoding an effector gene ranging from [0-1].

## `EffectorGeneP` parameter reference

```
python EffectorGeneP.py [-options] -g <genome fasta file> -t <transcripts in GFF3 format (gffread formatted)> -m <path to EffectorGeneP model files> -o <EffectorGeneP GFF3 output file>

Essential arguments are:
-g/--genome : genome FASTA file
-t/--transcript : assembled transcripts in GFF3 format (have to be gffread formatted)
-m/--model_files : path to EffectorGeneP model files, e.g. ./EffectorGeneP_Models/Fusarium_oxysporum_f_sp_lycopersici/
-o/--output : EffectorGeneP GFF3 output file

Options are:
-l/--length <int>: genes have to encode proteins with minimum length of <int> (default: 50 aas)
-s/--stringent : assume strand information of transcripts is always correct (not recommended for unstranded RNA-seq data)
-p/--padding <int> : add <int> bps to each transcript at both 5' and 3' ends to capture start/stops outside the transcript boundaries. Useful for low coverage data or overlapping UTRs (default: 200)
-c/--conservative : annotate genes conservatively

to enable sensitive search for genes encoding secreted proteins these two paths need to be provided:
--SIGNALP4 : path to signalp 4.1 executable
--TMHMM : path to tmhmm 2.0 executable

Advanced options are:
-d/--distance <int> : how far upstream to look for another translation start site, measured in encoded protein length (default: 70 aas, equates 210 nts)
-i/--intron <float> : for a multi-exon gene, average intron probability has to be >= <float> (default: 0.6, range from [0:1])

-h/--help : show brief help on version and usage
```


## FAQs

:question: Does `EffectorGeneP` give me a complete gene annotation?

Like TransDecoder (https://github.com/TransDecoder/TransDecoder), `EffectorGeneP` identifies candidate coding regions and genes within transcript sequences. Thus, it will currently only annotate genes that are expressed based on the transcripts you provide. Best practice for gene annotation is to merge evidence from multiple tools into a high-quality annotation. For targeted effector search such as in a defined genomic interval with high-quality transcriptional data from infection, `EffectorGeneP` can be used as a standalone tool. 

:question: Can I somehow provide protein homology evidence to `EffectorGeneP`?

Yes, it is possible to include homology evidence in the `EffectorGeneP` pipeline by combining transcripts with pseudo-transcripts derived from protein mappings to the reference genome. One could run a pipeline roughly like this:

```
# First, map reference proteins to genome
miniprot --outn=4 -G 3000 -J 10 --gff-only ${GENOME} {PROTEINS} > proteins.mini.gff3
# Only complete proteins with start, stop are kept 
gffread -J -g ${GENOME} proteins.mini.gff3 > proteins.mini.complete.gff3

# Merge with transcripts through gffcompare -C
gffcompare -C proteins.mini.complete.gff3 ${TRANSCRIPTS} -o evidence
gffread evidence.combined.gtf > combined.gff3
# Provide the file combined.gff3 file as input to `EffectorGeneP`
```

:question: `EffectorGeneP` is slow, how can I speed it up?

`EffectorGeneP` can be launched for each contig/scaffold/chromosome individually and all the individual output GFF3 files can then be concatenated into a full annotation GFF3 file. A simple job array will do the job and substantially speed up the annotation. For this, one could first use samtools to generate a list of contigs, split the input gff3 per contig and launch EffectorGeneP for each of these contig gff3 files:

```
# First, get the contig list for your genome
samtools faidx ${genome}
cat ${genome}.fai | cut -f1 > contigs.txt

# Now get the number of jobs that need to be run for a job array script
NUM_CONTIGS=$(wc -l < contigs.txt)
echo ${NUM_CONTIGS} "contigs will be annotated by EffectorGeneP."

# Then, launch your job array script like so:
sbatch --array=1-${NUM_CONTIGS}%50 EffectorGenes_JobArray.sh -a EffectorGeneP_input_transcripts.gff3 \
  -b ${OUTPUT_DIRECTORY} -c ${genome} -d ./EffectorGeneP_Models/Magnaporthe_oryzae/
```
You would then need to have a job array script (called EffectorGenes_JobArray.sh in this example) written up like so:

```
#!/usr/bin/env bash

#SBATCH --time=12:00:00
#SBATCH --job-name=EffectorGeneP_Array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5GB
#SBATCH --output=path_to_logfile_folder/effectorgenep_%A_%a.out

############################################################
module load python
############################################################
export OMP_NUM_THREADS=${SLURM_NTASKS}
helpFunction()
{
   echo ""
   echo "Usage: $0 -a INPUT_FILE -b OUTPUT_DIRECTORY -c GENOME_FASTA -d PATH_TO_MODELS"
   echo -e "\t-a the path to the INPUT GFF3 file"
   echo -e "\t-b the desired output directory for the EffectorGeneP annotations"      
   echo -e "\t-c the genome FASTA file"         
   echo -e "\t-d the path to the EffectorGeneP model files"         
   exit 1 # Exit script after printing help
}
############################################################
while getopts "a:b:c:d:" opt
do
   case "$opt" in
      a ) INPUT_FILE="$OPTARG" ;;
      b ) OUTPUT_DIRECTORY="$OPTARG" ;;
      c ) GENOME="$OPTARG" ;;
      d ) PATH_TO_MODELS="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done
############################################################
# Print helpFunction in case parameters are empty
if [ -z "$INPUT_FILE" ] || [ -z "$OUTPUT_DIRECTORY" ] || [ -z "$GENOME" ] || [ -z "$PATH_TO_MODELS" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi
############################################################
# Fetch the specific file path using the current array index line number
contig=$(sed -n "${SLURM_ARRAY_TASK_ID}p" contigs.txt)

awk -v search="${contig}" '$1 ~ /^#/ {print $0;next} {if ($1 == search) print}' $INPUT_FILE > $OUTPUT_DIRECTORY/"${contig}".gff3

echo ${contig}
mkdir OUTPUT_DIRECTORY

OUTPUT_FILE=${OUTPUT_DIRECTORY}"${contig}".EffectorGeneP.gff3
############################################################
# Begin script in case all parameters are correct
echo "-----------------"
echo "Use this input transcript assembly file (gffread-formatted):" $OUTPUT_DIRECTORY/"${contig}".gff3
echo "-----------------"
echo "Write the gene prediction files to this output file:" $OUTPUT_FILE
echo "-----------------"
echo "Use this genome FASTA file:" $GENOME
echo "-----------------"
echo "Use these EffectorGeneP model files:" $PATH_TO_MODELS
echo "-----------------"
############################################################
python EffectorGeneP.py -g ${GENOME} -t $OUTPUT_DIRECTORY/"${contig}".gff3 -m ${PATH_TO_MODELS} -o ${OUTPUT_FILE}
############################################################
```
After all jobs are finished, you would then need to simply concatenate the individual gff3 output files.

:question: There is no species model for my pathogen of interest, what do I do?

You can choose a closely related species to run `EffectorGeneP`. Have a look at the species models available here: https://effectorp.csiro.au/effectorgenep.html

:question: Can I apply `EffectorGeneP` to a bacterial pathogen genomes?

Not recommended at all as `EffectorGeneP` has only been trained on fungal genomes which have very different gene coding properties.


## Citation
If you use `EffectorGeneP` in your work please cite our preprint: https://www.biorxiv.org/content/10.64898/2026.05.04.722827v1

## License

Copyright (c) 2025 Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
This work is freely available for non-commercial scientific research, non-commercial education, or non-commercial research projects, under the CSIRO Non-Commercial License (https://github.com/JanaSperschneider/EffectorGeneP/blob/main/LICENCE). 

We welcome commercial enquiries and business partnership opportunities. Please contact: EnquiriesTeam@csiro.au

## Contact

Jana Sperschneider 

:office: https://people.csiro.au/s/j/jana-sperschneider 

:mailbox: jana.sperschneider@csiro.au


