#' Prepare data for `GAP`
#' @param ref_fsa Fasta file contains alignment of A and B genome ('-' stands for deletion).
#' @param ref_sam Sam file align B genome to A genome
#' @param alnA Sam file align reads to A genome
#' @param alnB Sam file align reads to B genome
#' @param ref_nameA Reference name of genome A which shows in file *refSam*.
#' With  ':' to separate genome name and base index and '-' to separate the begin (0-based)
#' and end index(begin + length of reference A, 1-based)
#' @param ref_nameB Reference name of genome B which shows in file *refSam*.
#' Similar way to define as genome A.
#' @param out_file File path and name of the outputted file
#' @param uni_geno_file File path and name of the universal fasta genome
#'
#' @note We could use MUMmer(and MUMmer3.20) to get sam file (ref.sam) and fasta file (ref.fa):
#' (sequence names in A.fa and B.fa should be named with the begin (0-based) and end (1-based)
#' index relative to the whole genomes)
#' nucmer -p ch1 --mum A.fa B.fa
#' delta-filter -1 ch1.delta > ch1_1to1.delta
#' MUMmer3.20/delta2maf ch1_1to1.delta > ch1_1to1.maf
#' cat ch1_1to1.maf | awk '/^s/{print ">" $2  "\\n" $7}' > ref.fa
#' nucmer --sam-long=refsam --mum A.fa B.fa
#' bioawk -c fastx '{  print "@SQ " "SN:" $name, "LN:" length($seq) }' < A.fa > header
#' cat header refsam > ref.sam

#' @useDynLib GAP r_make_aln
#' @import testthat
#' @export call_aln
#' @return Write the input datafile and universal reference genome file for `altraphase`
#' @examples
#' alignment <- system.file("extdata", "ref.fa", package = "GAP")
#' ref_sam <- system.file("extdata", "ref.sam", package = "GAP")
#' alnA <- system.file("extdata", "aln0A.sam", package = "GAP")
#' alnB <- system.file("extdata", "aln0B.sam", package = "GAP")
#' out_file <- system.file("extdata", "out.txt", package = "GAP")
#' uni_geno_file <- system.file("extdata", "uni.fa", package = "GAP")
#' call_aln(ref_nameA = "Genome_A:0-1000", ref_nameB = "Genome_B:0-1000",
#' ref_fsa = alignment, ref_sam = ref_sam, alnA = alnA, alnB = alnB,
#' out_file = out_file, uni_geno_file = uni_geno_file)

call_aln <- function(ref_nameA = NULL, ref_nameB = NULL, ref_fsa = NULL, ref_sam = NULL,
                     alnA = NULL, alnB = NULL, out_file = NULL, uni_geno_file = NULL) {
  checkmate::expect_file_exists(ref_fsa, access = "r")
  checkmate::expect_file_exists(ref_sam, access = "r")
  checkmate::expect_file_exists(alnA, access = "r")
  checkmate::expect_file_exists(alnB, access = "r")
  if (utils::tail(unlist(strsplit(ref_fsa, "[.]")), 1) != "fsa" &
      utils::tail(unlist(strsplit(ref_fsa, "[.]")), 1) != "fa" &
      utils::tail(unlist(strsplit(ref_fsa, "[.]")), 1) != "fasta")
    stop("The input ref_fsa has to be fasta file!")
  if (utils::tail(unlist(strsplit(ref_sam, "[.]")), 1) != "sam" &
      utils::tail(unlist(strsplit(alnA, "[.]")), 1) != "sam" &
      utils::tail(unlist(strsplit(alnB, "[.]")), 1) != "sam")
    stop("The input alnA, alnB, ref_sam must be sam file!")
  if (!dir.exists(sub('/[^/]*$', '', out_file)))
    stop("out_file path does not exist!")
  if (!dir.exists(sub('/[^/]*$', '', uni_geno_file)))
    stop("uni_geno_file path does not exist!")

  if (!is.loaded("r_make_aln", PACKAGE = "GAP"))
    dyn.load("../src/sync_data_r.so")

  .Call("r_make_aln", ref_nameA, ref_nameB, ref_fsa, ref_sam,
        alnA, alnB, out_file, uni_geno_file)
}

read_fasta <- function(datafile = NULL)
{
  checkmate::expect_file_exists(datafile, access = "r")
  if ((utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fasta") &&
      (utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fsa") &&
      (utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fa"))
    stop("The input datafile has to be fasta file!")

  if (!is.loaded("r_read_fasta", PACKAGE = "GAP")) {
    dyn.load("../src/sync_data_r.so")
  }
  res <- .Call("r_read_fasta", datafile)

  names(res) <- c("reads", "dim")
  if (length(res$reads) == res$dim[1] * res$dim[2])
    res$reads <- matrix(res$reads, ncol = res$dim[2], byrow = TRUE)

  return(res)
}

read_sam <- function(samfile = NULL, ref_name = NULL,
                     fastq_file = NULL, datafile = NULL) {
  checkmate::expect_file_exists(samfile, access = "r")
  if (utils::tail(unlist(strsplit(samfile, "[.]")), 1) != "sam")
    stop("The input datafile has to be sam file!")

  if (!is.loaded("./src/r_read_sam", PACKAGE = "GAP")) {
    dyn.load("../src/sync_data_r.so")
  }

  if(is.null(datafile) == TRUE | is.null(samfile) == TRUE)
    stop("Must input a data file and a samfile!")

  .Call("r_read_sam", samfile, ref_name, fastq_file, datafile)
}

read_fastq <- function(datafile = NULL)
{
  checkmate::expect_file_exists(datafile, access = "r")
  if (utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fastq")
    stop("The input datafile has to be fastq file!")

  if (!is.loaded("r_read_fastq", PACKAGE = "GAP")) {
    dyn.load("../src/sync_data_r.so")
  }

  res <- .Call("r_read_fastq", datafile)

  names(res) <- c("reads", "quality", "dim")
  res$reads <- matrix(res$reads, ncol = res$dim[2], byrow = TRUE)
  res$quality <- matrix(res$quality, ncol = res$dim[2], byrow = TRUE)

  return(res)
} #read_fastq
