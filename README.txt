MinION Pipeline 

This is a proposed pipeline to take long reads from MinION runs to then output assembled and aligned FASTA sequecenes based upon calling existing programmes.
The initial design will be based upon taking the 5kilobase Direct repeat region of Mycobacterium complex species and return consensus sequence (FASTQ and FASTA) and spoligotype.  

Rough pipeline: 
Input - FAST5 files; (future inputs would include: Genome assembly:Y/N; Length of seqeuence; reference sequence/genome)
Input > Guppy > demultiplexed basecalled FASTQ and FASTAs
FASTQ > SED or ORC > Barcoded sequences 
Sequences > Canu or Nanopolish > assemble sequences
assembled sequences > Nanopolish > Final Output 
(Output > SpoTyping-v2.1 > Spolgiotypes)

Will include multiple print statements along the way to indicated how the process is progressing (might include more print statements being added to source codes)