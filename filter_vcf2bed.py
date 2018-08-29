# AUTHOR: Johan Zicola
# DATE: 2018-04-09
# UPDATE:2018-04-11, add quality threshold to remove SNPs
# UPDATE:2018-04-12, Correct quality filter based on DP and GT. Works now.
# UPDATE:2018-05-04, add information in the output bed file  with REF, ALT, and quality field

# USAGE: Takes a VCF file as input (can be gzipped) and generates a BED file with the Ns
# in position that have REF=N, GT=. (missing data) OR/AND below a certain quality threshold defined
# By the user.
# To merge all contiguous Ns, just use bedtools merge tool (bedtools merge -i file.bed)
# I checked on a VCF file and N position in reference yields de facto a missing value in the call accession
# so the N reference should become N in sample w/o problems.
# Note that 1-based VCF positions are converted into 0-based BED positions


# Libraries needed. vcf should be downloaded on https://github.com/jamescasbon/PyVCF 
# Put folder in your bin and enter 'python setup.py test' to check if everything is running
# enter then 'python setup.py install --user'  to get the path to the library

import sys
import os
import argparse

# Check if PyVCF is installed or not
try:
    import vcf
except ImportError:
    raise ImportError("You must install PyVCF from https://github.com/jamescasbon/PyVCF in order to run this app")

    
def get_Ns(vcf_reader, sample, coverage, quality, exclude_deletions):

    # Check if only one sample in VCF file. Check size of list for record.genotype
    # if no argument 'sample' and only one sample in the VCF file, take it as 'sample'
    
    # Get list with samples from VCF file
    samples = vcf_reader.samples

    # Test if passed argument is empty or not
    if sample == None:
        if len(samples) > 1:
            sys.exit("The VCF file contains more than one sample, indicate sample name as second argument")
        else:
            sample = str(samples[0])
    
    # Check if sample provided is in the list
    try:
        test=samples.index(sample)
    except ValueError:
        sys.exit("The sample name provided is not in the VCF file. Provide valid sample name")


    # If no problem with sample, proceed with N call
    for record in vcf_reader:
        
        # Get chromosome
        chromosome = record.CHROM
        
        # Get position SNP. Note that VCF gives position in 1-based coordinates.
        # BED format uses 0-based coordinates and half-open (start-1, end)
        position = record.POS
        position_start = position-1
        position_end = position
        
               
        # Get info for the sample (if no optional argument, check if VCf file contains
        # only one sample. The output is given as a list, even when missing data (returns a .)
        # Note that althout GT returns a str, GQ and DP return lists in case they contain a value, and None 
        # in case they contain missing data (.)
        call = record.genotype(sample)
        
        # String giving the genotype of the sample
        GT = call['GT']
        
        # GQ and DP in Lists (unless = None), If the elements are in a list, extract them in a string
        GQ = call['GQ']
        if isinstance(GQ, (list,)):
            GQ = GQ[0]
        
        DP = call['DP']
        if isinstance(DP, (list,)):
            DP = DP[0]
        
        
        # Filter position
        
        # Note that a deletion with a quality below GQ and DP thresholds will
        # be considered as fitting the first loop statement, and therefore, will get
        # a coordinate of 1 bp. Since following nucleotide will also indicate the deletion
        # and be caught either in the first statement or the second.
        # In contrary, if the deletion has higher quality than GQ and DP threshold,
        # the full size of the deletion will be indicated in the bed file, unless
        # the option --exclude-deletions is on (in case deletions want to be incorporated in the
        # pseudogenome).
      
        # Call the nucleotide position in bed format if in these 3 cases:
        if call.data.GT == "." or GQ < quality or DP < coverage:
            # Replace the value 'None' by a dot for displaying missing data in output
            if GQ == None:
                GQ="."
            if DP == None:
                DP="."
            # ALT is always given as a list. Test if None and get a dot for output
            if record.ALT == [None]:
                record.ALT = "."
            # Print the output for the position matching the if conditions
            print str(chromosome)+"\t"+str(position_start)+"\t"+str(position_end)+"\t"+str(record.REF)+"\t"+str(record.ALT)+"\t"+str(GT)+":"+str(GQ)+":"+str(DP)
        
         
        # If the option exclude_deletions is not selected, and the GT as a value so that
        # an alternative allele (GT > 0) is there and has DP > coverage and GQ > quality 
        # (did not match the previous if statement. Call the reference allele of 
        # the sample and test whether it is smaller than the reference allele 
        # (hence there is a deletion)
        # if yes, get the position + the length of the deletion and include as bed entry

        elif not exclude_deletions and GT > "0" :
            index = int(call.data.GT)
            # Get the sequence of the alternative allele for the sample
            ALT = record.ALT[index-1]
            # Check if there is something in ALT and if string ALT is smaller than REF
            #if (ALT != None and len(ALT) < len(record.REF)) and (GQ < quality or DP < coverage):
            if ALT != None and len(ALT) < len(record.REF):
                length_deletion = len(record.REF)
                position_end = position_start + length_deletion
                print str(chromosome)+"\t"+str(position_start)+"\t"+str(position_end)+"\t"+str(record.REF)+"\t"+str(ALT)+"\t"+str(GT)+":"+str(GQ)+":"+str(DP)


                
# Function dealing with the arguments. More convenient to separate in its own block
def handle_argument():

    #Create parser
    parser = argparse.ArgumentParser(description='AUTHOR: Johan Zicola\nDATE: 2018-04-11\nTakes a VCF file as input (can be gzipped)\
            and generates a BED file with the Ns in position that have REF=N, GT=. (missing data)\
            OR/AND below a certain quality threshold defined by the user. To merge all contiguous\
            Ns, just use bedtools merge tool (bedtools merge -i file.bed).')

    #name of argument needed, the name after -- is used for call from args object (echo, not e)
    parser.add_argument("-i","--input_vcf", help="VCF file", type=str) 
    # Note that I can use underscore in argument name but no hyphens
    parser.add_argument("-s","--sample", help="Sample name", default=None, type=str)
    parser.add_argument("-c","--coverage", help="Threshold coverage (DP) for the SNPs (default=3)", default=3, type=int)
    parser.add_argument("-q","--quality", help="Threshold quality (GT) for the SNPs (default=25)", default=25, type=int)
    parser.add_argument("-e","--exclude-deletions", help="Exclude deletions from filtering if -e added", action="store_true")

    # Display help message if no argument is provided
    if len(sys.argv)==1:
            parser.print_help(sys.stderr)
            sys.exit(1)
    
    args=parser.parse_args()
    
    return args


   
def main():
    
    args = handle_argument()
    # The PyVCF package handles compression of the file (if bgzipped)
    vcf_reader = vcf.Reader(open(args.input_vcf, 'r'))
    
    #get_all_nucleotides()
    get_Ns(vcf_reader, args.sample, args.coverage, args.quality, args.exclude_deletions)
    
        
if __name__ == "__main__":
    sys.exit(main())


