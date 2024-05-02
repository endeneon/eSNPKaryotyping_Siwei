# eSNPKaryotyping_Siwei
Rewrittened eSNPKaryotyping script using GATK 4.1.8.1 (3.x are not supported any more);
Refer to https://github.com/BenvenLab/eSNPKaryotyping for the original scripts and citations;

Required but not provided (for the size)
- R >= 4.2.1;
- GATK [4.1.8.0, 4.3.0.0). Currently not support GATK > 4.3.0.0 due to syntax changes;
- Human GRCh38 genome fasta file, such as the one provided by the 10x Genomics (refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa);
- Humab db150 common SNP files separated by chromosomes, available at ftp.ncbi.nlm.nih.gov/snp. Check if the dbSNP database shares the same chromosome nomination rules (e.g.: chr1 but not 1) as the genome reference fasta and are both in GRCh38;
- After git clone, some environment-specific paths need to be revised to make accordance with the local environment.
