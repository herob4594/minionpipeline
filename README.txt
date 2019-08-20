MinION Pipeline 

Requirements:
Guppy
Canu
Nanopolish
Minimap2
Samtools
Blast
Python(2.7)
Blast

This is a proposed pipeline to take long reads from MinION runs to then output assembled and aligned FASTA sequecenes based upon calling existing programmes.
The initial design will be based upon taking the 5kilobase Direct repeat region of Mycobacterium complex species and return consensus sequence (FASTQ and FASTA) and spoligotype.  

Rough pipeline: 
Raw FAST5 (input1) > Guppy_basecaller > FASTQ > canu > draft.assembly
draft assembly, FASTQ, FAST5 > nanopolish > minimap2 > samtools > nanopolish > polished /consensus sequence 
consensus seqeuence > SpoTyping-v2.1 > Spoligotype 

Will include multiple print statements along the way to indicated how the process is progressing (might include more print statements being added to source codes)