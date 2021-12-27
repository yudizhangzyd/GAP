`GAP` is an R package developed to simultaneously genotyping and phasing alleotetroploid. The package can be installed from GitHub using devtools and then loaded in the usual way.

Installation
------------

Our package calls function from [mnlogit](https://www.rdocumentation.org/packages/mnlogit/versions/1.2.6/topics/mnlogit), please install it first.

`GAP` can be installed from the GitHub repository using the devtools package.

``` r
devtools::install_github("yudizhangzyd/GAP")
```

Once installed, the package can be loaded as usual.

``` r
library(GAP)
```

Functions
---------

The package contains the following functions.

## call_aln

`call_aln` prepares data for our main function `altraphase`.

#### Input (required):
		
1. ref_sam
		
SAM file containing alignment of references (selected regions).

   - The names of the selected regions needs to follow: The names of regions need to contain start and end position relative to the whole genome and the genome name. Using ':' to seperate original genome name and region index, '-' to seperate start and end position. For example, chr1:0-10 means we will genotype and phase from position 1 to 10 (1 based) in chr1 genome.

2. ref_fsa
	
Reference fasta files contain BOTH aligned selected references (names of references must match the names in ref_sam). 

  -  The indels in alignments must be represented as `-`.
  Users can use functions `nucmer`, `delta2maf`, `nucmer` from [MUMmer4](https://github.com/mummer4/mummer) and [MUMmer3.20](https://sourceforge.net/projects/mummer/files/mummer/3.20/) to obtain both files for ref_sam and ref_fsa. Example below:
  ```
  nucmer -p ref --mum A_ref.fa B_ref.fa 
  MUMmer3.20/delta2maf ref.delta > ref.maf
  awk '/^s/{print ">" $2 "\n" $7}' ref.maf > ref.fa
  nucmer --sam-long=ref --mum A_ref.fa B_ref.fa
  ```
  This will give ref.fa and ref.sam required for our package.
  
Note: for `delta2maf` function in MUMmer3.20, you can install it as:
  ```
  svn checkout svn://svn.code.sf.net/p/mugsy/code/trunk/MUMmer3.20/
  cd MUMmer3.20/src/tigr
  make delta2maf
  ```
  
3. alnA and alnB 

Alignment of reads to selected subgenomic regions (A_ref.fa and B_ref.fa in previous example).
   
#### Output: 

out_file and uni_geno_file, contain the information of the alignments of reads and the alignment of references respectively.

#### Example: 
```
alignment <- system.file("extdata", "ref.fa", package = "GAP")
ref_sam <- system.file("extdata", "ref.sam", package = "GAP")
alnA <- system.file("extdata", "aln0A.sam", package = "GAP")
alnB <- system.file("extdata", "aln0B.sam", package = "GAP")
out_file <- system.file("extdata", "out.txt", package = "GAP")
uni_geno_file <- system.file("extdata", "uni.fa", package = "GAP")
call_aln(ref_nameA = "Genome_A:0-1000", ref_nameB = "Genome_B:0-1000", ref_fsa = alignment, ref_sam = ref_sam, alnA = alnA, alnB = alnB, out_file = out_file, uni_geno_file = uni_geno_file)
```
   
## altraphase

`altraphase` main function for genotyping and phasing alleotetraploid individual. It takes the output of `call_aln` and returns a list contain the phased haplotypes, SNPs and SNP locations and other auxiliary values.

#### Example:

Users can set parameters such as `max_iter` (max iteration, default 20) and `tol` (onvergence tolerance, default 1e-05)
```
datFile <- system.file("extdata", "out.txt", package = "GAP")
refFile <- system.file("extdata", "uni.fa", package = "GAP")
res <- altraphase(datafile = datFile, alignment = refFile, max_iter = 10)
```


