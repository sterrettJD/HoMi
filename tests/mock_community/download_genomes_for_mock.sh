# E. coli K-12 substr. MG1655
mkdir -p e_coli
cd e_coli
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000005845.2/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000005845.2.zip" -H "Accept: application/zip"
unzip GCF_000005845.2
cd ..


# Clostridium beijerinckii str. NCIMB 8052
mkdir -p c_beijerinckii
cd c_beijerinckii
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000016965.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000016965.1.zip" -H "Accept: application/zip"
unzip GCF_000016965.1
cd ..


# Faecalibacterium prausnitzii M21/2
mkdir -p f_prausnitzii
cd f_prausnitzii
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000154385.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000154385.1.zip" -H "Accept: application/zip"
unzip GCF_000154385.1
cd ..


# Human pangenome
mkdir -p human_pangenome
cd human_pangenome
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph/CHM13v11Y.fa.gz
gunzip CHM13v11Y.fa.gz
cd ..