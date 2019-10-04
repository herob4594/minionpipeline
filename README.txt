README.txt
MinION Pipeline 

Components:
pipeline_script2.sh (bash script)


Requirements:

Guppy
canu 1.8
Nanopolish
Minimap2 (minimap2-2.17_x64-linux)
Samtools
Blast
Python2.7
Spotyping-2.1 (uses Blast)


IMPORTANT!
All these programmes have been assumed to be within the same directory (see code for more details), and the user will be asked to state the directory all the above programmes are stored in. 

This is a proposed pipeline to take long reads from MinION runs to then output assembled and aligned FASTA sequences based upon calling existing programmes.
The initial design will be based upon taking the 5kilobase Direct repeat region of Mycobacterium complex species and return a consensus sequence (FASTQ and FASTA) and spoligotype.  
It includes multiple print statements along the way to indicate how the process is progressing. 

Rough pipeline: 
MinION run with PCR amplified 5kb direct repeat region from M. bovis organism
Raw FAST5 (input1) > Guppy_basecaller > FASTQ > canu > draft.assembly
draft assembly, FASTQ, FAST5 > nanopolish > minimap2 > samtools > nanopolish > polished /consensus sequence 
consensus sequence > SpoTyping-v2.1 > Spoligotype 

All steps included:
This is included as it is recommended to try each of these steps individually if the script is not working, as well as giving insight into each step and the pipeline.

PATH_TO_FAST5/ and SUFFIX will be user defined inputs

The following also acts as a guide on how to produce a consensus sequence (and then a spoligotype) from raw MinION data.

guppy_basecaller -r --input_path /PATH_TO_FAST5/ --save_path ./guppy_SUFFIX --config dna_r9.4.1_450bps_hac.cfg --qscore_filtering --min_qscore 7 --cpu_threads_per_caller 8 --num_callers 1
#basecalling

/opt/canu-1.8/Linux-amd64/bin/canu -p draft_assembly_SUFFIX -d ./canu genomeSize=5k -nanopore-raw ./guppy_SUFFIX/pass/fastq
#canu producing a draft assembly of the 5 kilobase direct repeat region

/opt/nanopolish/nanopolish index -d PATH_TO_FAST5 ./guppy_SUFFIX/pass/fastq
#index fast5 reads with the fastq reads 

/opt/minimap2-2.17_x64-linux/minimap2 -ax map-ont -t 8 ./canu/draft_assembly_SUFFIX.contigs.fasta ./guppy_SUFFIX/pass/fastq > alignment_SUFFIX.sam 
#minimap2 alignment to sam file

samtools view -S -b alignment_SUFFIX.sam > alignment_SUFFIX.bam
#sam to bam file conversion (samtools sorting need bam file input)

samtools sort -mof alignment_SUFFIX.bam alignment_SUFFIX.sorted
#sorting the bam file with samtools 

samtools index alignment_SUFFIX.sorted.bam
#samtools bai index file made

/opt/nanopolish/nanopolish variants --consensus -o polished_SUFFIX.vcf -r ./guppy_SUFFIX/pass/fastq -b alignment_SUFFIX.sorted.bam -g ./canu/draft_assembly_SUFFIX.contigs.fasta
#Vcf polished consensus sequence 

/opt/nanopolish/nanopolish vcf2fasta --skip-checks -g ./canu/draft_assembly_SUFFIX.contigs.fasta polished_SUFFIX.vcf > polished_SUFFIX.fasta
#Now consensus in fasta format

python2.7 /opt/SpoTyping-2.1/SpoTyping-v2.1-commandLine/SpoTyping.py --seq polished_SUFFIX.fasta -o SUFFIX.spoligotype
#Spoligotype then taken from the consensus direct repeat region 


A paper displaying the full methods of the whole process (from sampling of bacteria to returned in silico spoligotype) will be released soon. 
