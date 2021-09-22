# Genome assembly pipeline for *Drosophila pseudoobscura* (ST) and *Drosophila persimilis* (M40) - Extended methods.
This is a genome assembly pipeline used (from paper) for the *D. pseudoobscura* (ST) and *D. persimilis* (M40).
This pipeline implements a hybrid assembly using Illumina short reads and PacBio CLR reads.

## Step 1
Generation of a PacBio-only, gap-filled and polished assembly: this assembly was construted directly by Pacific Biosciences using HGAP-Arrow with default parameters (please refer to the [HGAP github page](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/HGAP) for details). The assembly statistics for each species are available at the Supplementary_Figures_Tables.docx file (Table S13).

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
