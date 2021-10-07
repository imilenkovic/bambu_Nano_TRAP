library(bambu)

### Load .bam files of Input samples placed in the folder in which the script is being run

bamFiles <- Rsamtools::BamFileList(c("Heart_WT_adult_rep1.bam", "Heart_WT_adult_rep2.bam", "Heart_WT_adult_rep3.bam", "Heart_KO_adult_rep1.bam", "Heart_KO_adult_rep2.bam", "Heart_KO_adult_rep3.bam"))

### Load the gtf annotation (in order to match chromosome names in .bam files we had to rename "Chr 1"  to "1" in the gencode annotation; that's why the file is called "no_chr")

gtf.file <- "no_chr_gencode.vM25.basic.annotation.gtf"
bambuAnnotations <- prepareAnnotations(gtf.file)

### Load the reference - all non-protein coding sequences were masked from the GRC38 assembly

### Masking the reference was done as follows:

#### remove all protein coding genes to keep only the ones you want to mask
grep -v "protein_coding" gencode.vM25.basic.annotation.gtf > masking.annotation.gtf

#### convert gtf to bed - THIS DIDNT WORK, but it's still fine if it's gtf
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' masking.annotation.gtf | gtf2bed - > masking.annotation.bed

#### convert "chr1" to "1" (because it's like this in the genome fasta file)
sed 's/chr//g' masking.annotation.bed > ivanito.bed

#### mask using the created bed
bedtools maskfasta -fi Mus_musculus.GRCm38.dna.primary_assembly.fa -fo genome_masked.fa -bed ivanito.bed

### Load the reference

genomeSequence <- "genome_masked.fa"

### Run bambu

se <- bambu(reads = bamFiles, annotations = bambuAnnotations, genome = genomeSequence,
    ncore = 4, discovery = FALSE)

colData(se)$condition <- as.factor(c("WT", "WT", "WT", "KO", "KO", "KO"))
seGene <- transcriptToGeneExpression(se)

#save.dir <- "/nfs/users2/enovoa/imilenkovic/software/bambu/"

### Run DESeq on the bambu object

dds.deseq <- DESeq(dds)
deGeneRes <- DESeq2::results(dds.deseq, independentFiltering = FALSE)
head(deGeneRes[order(deGeneRes$padj), ])

summary(deGeneRes)

### WRITE TO RESULTS

write.csv(deGeneRes, "Bambu_no_novel_prediction_DESeq2_IP_24082021.csv", row.names = FALSE)
