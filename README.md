# Genome assembly pipeline for *Drosophila pseudoobscura* (ST) and *Drosophila persimilis* (M40) - Extended methods.
This section shows the genome assembly pipeline used (from paper) for the *D. pseudoobscura* (ST) and *D. persimilis* (M40).
This pipeline implements a hybrid assembly using Illumina short reads and PacBio CLR reads.

## Step 1
Generation of a PacBio-only, gap-filled and polished assembly: this assembly was construted directly by Pacific Biosciences using HGAP-Arrow with default parameters (please refer to the [HGAP github page](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/HGAP) for details). The assembly statistics for each species are available at the Supplementary_Tables.xlsx file (Table S13).

### Step 1.1
- Gap-filling of the PacBio-only assemblies using [PbJelly](https://sourceforge.net/projects/pb-jelly/files/).

Commands:

>`Jelly.py mapping jellyProtocol.xml`

>`Jelly.py support jellyProtocol.xml`

>`Jelly.py extraction jellyProtocol.xml`

>`Jelly.py assembly jellyProtocol.xml -x --nproc=8`

> NOTE: Please refere to the jellyProtocol.xml for full parameter details.

### Step 1.2
- Polishin using [Pilon](https://github.com/broadinstitute/pilon/releases/tag/v1.22)

Commands:

> `bwa mem -t 5 st_pbjelly.fasta file_1.fastq file_2.fastq > pbjelly.sam 2> stderror.txt`

> `samtools view -bS pbjelly.sam > pbjelly.bam 2> stderror_samtoolsview.txt`

> `samtools sort -@ 5 -o pbjelly_sorted.bam pbjelly.bam 2> stderror_samtoolssort.txt`

> `samtools view -b -F 12 pbjelly_sorted.bam > pbjelly_mapped_sorted.bam 2> stderror_samtoolsmap.txt`

> `samtools index pbjelly_mapped_sorted.bam 2> stderror_samtoolsindex.txt`

> `java -Xmx30G -jar pilon-1.22.jar --genome pbjelly.fasta --frags pbjelly_mapped_sorted.bam --threads 5 --changes 2> stderror_pilon.txt`

## Step 2.
Generation of a hybrid assembly using CLR reads and Illumina paired-end reads using [DBG2OLC](https://github.com/yechengxi/DBG2OLC).

Commands:

> `SparseAssembler LD 0 k 51 g 15 NodeCovTh 1 EdgeCovTh 0 GS 171281433 i1 file_1.fastq o1 file_2.fastq > sparseAssembler_test1.log`

> `SparseAssembler LD 1 k 51 g 15 NodeCovTh 2 EdgeCovTh 1 GS 171281433 i1 st_1.fastq o1 st_2.fastq > sparseAssembler_test1.log`

> `DBG2OLC k 17 KmerCovTh 2 MinOverlap 20 AdaptiveTh 0.002 LD1 0 MinLen 200 Contigs Contigs.txt RemoveChimera 1 f pacbio.fasta`

## Step 3.
Final scaffolding step using [quickmerge](https://github.com/mahulchak/quickmerge)

Commands:

>`nucmer -l 100 -p out -t 10 pacbiopolished.fasta hybrid.fasta 2> stderror_nucmer.txt`

> `delta-filter -i 95 -r -q out.delta > out.rq.delta 2> stderror_deltafilter.txt`

> `quickmerge -d out.rq.delta -q hybrid.fasta -r pacbiopolished.fasta -hco 5.0 -c 1.5 -l n -ml m 2> stderror_quickmerge.txt`

# Structural variation analysis.
This section shows the followed approach for the structural variation analysis that includes INDEL and CNV calling using CLR reads and whole genome alignments. 


