Make pseudogenome
=====

# Usage
Generate pseudogenomes derived from a reference fasta file and a VCF file containing one or several samples. The output file is in the fasta format. 

Note that the python script filter_vcf2bed.py was originally written to process VCF file from homozygous species (*i.e. Arabidopsis thaliana* in our lab). Some modifications for phased heterozygous calls may be required.

# Softwares/scripts required
* bcftools v1.2
* vcftools v0.1.14
* python script filter_vcf2bed.py (written for Python2.7)
* tabix

# Pipeline
Two methods are suggested here to generate a pseudogenome. The first method is quicker and replaces all positions of low quality or missing data in the VCF file by the reference allele. The second method involves the use of a mask which allows to replace unwanted position (e.g. low quality call, missing data, deletions, ...) by Ns. The latter method yields a more correct pseudogenome *sensu stricto* but generates pseudogenomes with often lower mappability as the assumption of the reference allele is often correct. For  applications related to mapping short reads (BS-seq, DNase-seq, ChIP-seq, ...), the first method is advised.

## Method 1 - without mask
### Step 1: Subset VCF file with indels
This step allows to process faster the downstream steps as computational time increases with the size of the VCF file used. Only positions containing at least 1 alternative allele is kept (--min-alleles 2) allowing to generate a much smaller VCF file. The tool *bcftools consensus* assumes reference allele when no data are available so that positions removed will be *de facto* considered as reference. The vcf file should be compressed with bgzip and indexed with tabix. The samples_id.txt should contain the accessions ID (matching in VCF), one line per accession. If needed to check samples ID, type "bcftools query -l <file.vcf.gz>". One can directly  filter the VCF file to keep only positions with a certain quality threshold. If the chromosome lengths should stay identical to the reference in the pseudogenome, indels should be removed from the VCF file.  

Subsetting by keeping positions with GQ >= 30 and DP >= 5 and removing indels:
```
vcftools --gzvcf <subset_vcf_file.gz> --remove-indels --minGQ 30 --minDP 5  --recode --recode-INFO-all --out subset_vcf_file_wo_indels 
```
Subsetting by keeping positions with GQ >= 30 and DP >= 5 and keeping indels:
```
vcftools --gzvcf <subset_vcf_file.gz> --minGQ 30 --minDP 5  --recode --recode-INFO-all --out subset_vcf_file_with_indels 
```

gzip and create an index:
```
bgzip subset_vcf_file_wo_indels.vcf
tabix subset_vcf_file_wo_indels.vcf.gz
```

### Step 2: Generate pseudogenome
Before performing this command, verify that the chromosome nomenclature in the reference fasta file is the same than in the VCF file:
```
# See chromosome names in the VCF file
bcftools view -h subset_vcf_file_wo_indels.vcf.gz | grep "##contig" -

# See chromosome names in reference fasta file
grep ">" reference_fasta.fa
```

Generate the pseudogenome:
```
bcftools consensus <subset_vcf_file_wo_indels.vcf.gz> --sample sample_name --fasta-ref <reference_fasta.fa> > pseudogenome_sample_name.fa
```
Use the command in a while loop to generate multiple pseudogenomes. For example, if all samples contained in the VCF file needs to be converted into a fasta pseudogenome:
```
while read i; do
    bcftools consensus <subset_vcf_file_wo_indels.vcf.gz> --sample $i --fasta-ref <reference_fasta.fa> > pseudogenome_${i}.fa"
done <<< $(bcftools query -l <subset_vcf_file_wo_indels.vcf.gz>)
```

## Method 2 - with mask
### Step 1: Subset VCF file with indels
This step allows to process faster the downstream steps as computational time increases with the size of the VCF file used. Only positions containing at least 1 alternative allele is kept (--min-alleles 2) allowing to generate a much smaller VCF file. The tool *bcftools consensus* assumes reference allele when no data are available so that positions removed will be *de facto* considered as reference. The vcf file should be compressed with bgzip and indexed with tabix. The samples_id.txt should contain the accessions ID (matching in VCF), one line per accession. If needed to check samples ID, type "bcftools query -l <file.vcf.gz>

```
vcftools --keep <samples_id.txt> --gzvcf <file.vcf.gz> --min-alleles 2 --recode --out subset_vcf_file
```
gzip and create an index:

```
bgzip subset_vcf_file.vcf
tabix subset_vcf_file.vcf.gz
```

### Step 2: Generate a VCF subset of your accessions of interest without indels
This allows to obtain a pseudogenome with the same size than the reference accessions. Note that if the difference in size of the chromosome is not an issue, this step can be skipped.

Create a new vcf file without indels
```
vcftools --gzvcf <subset_vcf_file.gz> --remove-indels --recode --recode-INFO-all --out subset_vcf_file_wo_indels
```

gzip and create an index:
```
bgzip subset_vcf_file_wo_indels.vcf
tabix subset_vcf_file_wo_indels.vcf.gz
```

Note that the deletions removed from the VCF file will be replaced by reference  allele when using bcftools consensus tool. Since no reads should be mapping at the deletions (if resequencing is performed), this does not create bias in alignment apart from having reads being multireads in the reference but uniquely mapped in the pseudogenome (due to deletions). 

### Step 3: Generate mask

The mask is a bed file which contains all regions to be ignored when the pseudogenome is made (nucleotides are replaced by Ns). It contains positions with low quality score according to user's defined thresholds. Deletions can be added to the mask if needed to be replaced by Ns. This enable to avoid masking these deletions and therefore integrate them in the pseudogenome. Note that this option will generate pseudogenomes of unequal sizes.

If the pseudogenome must have the same chromosome size than the reference, the deletions should be replaced by Ns (default mode) and the insertions should be ignored :
```
python filter_vcf2bed.py -s sample_name -i <subset_vcf_file.vcf.gz> > mask_sample.bed
```

*NB: Insertions are usually absent in SNP VCF files as common SNP caller do not deal well with insertions or structural variants in general.*

The default quality parameters are GQ >= 25 and DP >= 3. If for example GQ >= 30 and DP >= 5 are required:
```
python filter_vcf2bed.py -s sample_name -i <subset_vcf_file.vcf.gz> -q 30 -c 5 > mask_sample.bed
```

In the case deletions should be kept if passing the quality test (resulting in unequal genome size), the option *--exclude-deletions* should be added to the command so they are not added to the mask bed file:
```
python filter_vcf2bed.py --exclude-deletions -s sample_name -i <subset_vcf_file.vcf.gz> > mask_sample_with_deletions.bed
```

## Step 4: merge bad quality call and deletions in one bed file to use as mask
This step is optional but yiels much smaller bed file as contiguous positions are merged.

Merge mask_sample.bed (or mask_sample_with_deletions.bed) file:
```
bedtools merge -i mask_sample.bed > mask_sample.merged.bed
```

Rename and overwrite initial unmerged bed file:
```
mv mask_sample.merged.bed mask_sample.bed
```

Use 'rename -f 's/\.bed.merged.bed$/.bed/' *merged.bed' if several bed files need to be renamed

### Step 5: Generate pseudogenome
Two alternatives are possible:

#### 1. Generate pseudogenomes without deletions
The mask will enable to replace by Ns all bases with a bad quality call and mask the deletions.
```
bcftools consensus <subset_vcf_file_wo_indels.vcf.gz> --sample sample_name --fasta-ref <reference_fasta.fa> --mask <mask_sample.bed> > pseudogenome_sample_name.fa
```

### 2. Generate pseudogenome with deletions
In this case, the mask bed file was generated using the option -e (exclude deletions) to allow deletions passing a certain quality to be included in the pseudogenome:
```
bcftools consensus <subset_vcf_file.vcf.gz> --sample sample_name --fasta-ref <reference_fasta.fa> --mask <mask_sample_with_deletions.bed> > pseudogenome_with_deletions_sample_name.fa
```

# Examples

Check [examples](https://github.com/johanzi/make_pseudogenome/blob/master/examples_filter_vcf2bed/script.sh) to understand how "bcftools consensus" and filter_vcf2bed.py are working.


## Authors

* **Johan Zicola** - [johanzi](https://github.com/johanzi)

Script written on 2018-05-29

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details


