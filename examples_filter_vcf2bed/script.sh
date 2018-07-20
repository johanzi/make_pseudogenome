# Exemples of how to use bcftools consensus tool and the filter_vcf2bed.py script
# to generate pseudogenomes

# VCF file. The fields indicated are required. In addition, bcftools consensus
# takes as input only vcf compressed with bgzip and indexed with tabix

cat test.vcf
##fileformat=VCFv4.0
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FILTER=<ID=q25,Description="Quality below 25">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##contig=<ID=chloroplast>
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	13582	6911
chloroplast	1	.	A	C,TA	14	q25	.	GT:GQ:DP	1:80:82	2:30:33
chloroplast	2	.	T	G,GG	17	q25	.	GT:GQ:DP	1:14:97	2:30:1
chloroplast	3	.	G	A	17	q25	.	GT:GQ:DP	1:34:97	0:30:20
chloroplast	4	.	A	AAAAAAA	24	q25	.	GT:GQ:DP	1:26:70	1:30:5
chloroplast	5	.	C	T	40	q25	.	GT:GQ:DP	.:.:.	.:.:.
chloroplast	6	.	C	T	40	q25	.	GT:GQ:DP	.:22:2315	0:29:8
chloroplast	7	.	C	G	40	q25	.	GT:GQ:DP	.:28:2315	0:29:8

# Compress vcf file
bgzip test.vcf -c > test.vcf.gz

# Index compressed vcf file
tabix test.vcf.gz

# Reference fasta file seq (require fasta header)
cat test.fa
>chloroplast
ATGACCC

# Make consensus w/o mask for the sample '13582'
bcftools consensus test.vcf.gz --sample 13582 --fasta-ref test.fa

>chloroplast
CGAAAAAAAACCC


# Make a mask
python ~/SCRIPTS/filter_vcf2bed.py -i test.vcf.gz -s 13582 > mask.bed

cat mask.bed
chloroplast     1       2       T       [G, GG] 1:14:97
chloroplast     4       5       C       [T]     .:.:.
chloroplast     5       6       C       [T]     .:22:2315
chloroplast     6       7       C       [G]     .:28:2315


# Remake consensus with a mask
bcftools consensus test.vcf.gz --sample 13582 --fasta-ref test.fa --mask mask.bed

>chloroplast
CNAAAAAAAANNN


# Change quality threshold for mask
python ~/SCRIPTS/filter_vcf2bed.py -i test.vcf.gz -q 30 -s 13582 > mask.bed
 
chloroplast     1       2       T       [G, GG] 1:14:97
chloroplast     3       4       A       [AAAAAAA]       1:26:70
chloroplast     4       5       C       [T]     .:.:.
chloroplast     5       6       C       [T]     .:22:2315
chloroplast     6       7       C       [G]     .:28:2315

bcftools consensus test.vcf.gz --sample 13582 --fasta-ref test.fa --mask mask.bed

# The insertion is know replaced by 1 N
>chloroplast
CNANNNN
