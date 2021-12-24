`GAP` is an R package developed to simultaneously genotyping and phasing alleotetroploid. The package can be installed from GitHub using devtools and then loaded in the usual way.

Installation
------------

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

#### call_aln

`call_aln` Prepare data for `altraphase`.

Input (required):
		
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
  This will give ref.fa and ref.sam reequired for our package.
  
Note: for `delta2maf` function in MUMmer3.20, you can install it as:
  ```
  svn checkout svn://svn.code.sf.net/p/mugsy/code/trunk/MUMmer3.20/
  cd MUMmer3.20/src/tigr
  make delta2maf
  ```
  
3. alnA and alnB 
   Alignment of reads to selected subgenomic regions (A_ref.fa B_ref.fa in previous example).
#### altraphase

`altraphase` main function for genotyping and phasing alleotetraploid individual.




