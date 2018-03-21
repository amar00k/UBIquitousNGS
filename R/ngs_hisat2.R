#' NGS mapping with Hisat2
#'
#' @import UBIquitous
#'
#' @param reads1.pathnames samples to map
#' @param sample.ids basenames of ouput files
#' @param hisat2.command command to run hisat2
#' @param hisat2.build.command command to run hisat2-build
#' @param genome.fasta genome fasta reference
#' @param genome.index.name genome hisat2 index basename
#' @param output.dir directory to store output files
#' @param MAX.THREADS maximum number of threads to use
#' @param reads2.pathnames
#'
#' @return
#' @export
#'
#' @examples
ngs_hisat2 <- function(reads1.pathnames,
                       sample.ids,
                       hisat2.command,
                       hisat2.build.command,
                       genome.fasta,
                       genome.index.name,
                       output.dir,
                       MAX.THREADS=8,
                       reads2.pathnames = rep(NA, length(reads1.pathnames)),
                       title="") {
  MODULE.DESCRIPTION <- "This step maps the sequenced reads to a reference genome
  using [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)."

  # reads2.pathnames <- ifelse(is.null(reads2.pathnames), rep(NA, length(reads1.pathnames)), reads2.pathnames)

  # build the index
  cmd.build <- paste(hisat2.build.command,
                     genome.fasta,
                     genome.index.name)

  file.check <- paste0(genome.index.name, ".1.ht2")

  # check for index and build
  res <- integer(0)
  if (!file.exists(file.check)) {
    message <- system(paste(cmd.build, "2>&1"), intern=TRUE)
    res <- attr(message, "status")
  }

  if (length(res) != 0) {
    cat('<div class="alert alert-danger">')
    cat('<b> Error building genome index: </b><br>')
    cat(message, sep="<br>")
  }

  # generate mapping commands
  sam.to.bam.commands <- function(target.dir, sample.id, deleteSam=TRUE) {
    sam <- file.path(target.dir, paste0(sample.id, ".sam"))
    bam <- file.path(target.dir, paste0(sample.id, ".bam"))
    bam.sorted.basename <- file.path(target.dir, paste0(sample.id, ".sorted"))
    bam.sorted <- file.path(target.dir, paste0(sample.id, ".sorted.bam"))

    c(paste("samtools", "view", "-bS", sam, ">", bam),
      paste("samtools", "sort", bam, bam.sorted.basename),
      paste("samtools", "index", bam.sorted),
      paste("rm", bam, ifelse(deleteSam, sam, "")))
  }

  make.hisat2.commands <- function(reads1, reads2, sample.id, target.dir, deleteSam=TRUE) {
    flag.type <- "-q"

    flag.reads <- ifelse(is.na(reads2),
                         paste("-U", reads1),
                         paste("-1", reads1, "-2", reads2))

    map.cmd <- paste(hisat2.command,
                     flag.type,
                     # "-p", MAX.THREADS,
                     "--new-summary --summary-file", file.path(target.dir, paste0(sample.id, ".log")),
                     "-x", genome.index.name,
                     flag.reads,
                     ">", file.path(target.dir, paste0(sample.id, ".sam")))

    paste(c(map.cmd, sam.to.bam.commands(target.dir, sample.id)), collapse="\n")
  }

  target.files <- file.path(output.dir, paste0(sample.ids, ".sorted.bam"))

  map.cmds <- unlist(lapply(1:length(sample.ids),
                            function(i) make.hisat2.commands(reads1.pathnames[i],
                                                             reads2.pathnames[i],
                                                             sample.ids[i],
                                                             output.dir)))
  # run the mapping
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

  w <- which(!file.exists(target.files))
  map.res <- unlist(parallel::mclapply(map.cmds[w], system, mc.cores = MAX.THREADS))

  # multiqc on alignments
  multiqc.file <- file.path(output.dir, "multiqc_report.html")

  if (!file.exists(multiqc.file)) {
    system(paste("multiqc", "-o", output.dir, output.dir))
  } else {
    # check that all samples are processed by multiqc
    mqc <- read.table(file.path(output.dir, "multiqc_data", "multiqc_hisat2.txt"), sep="\t", header=TRUE)

    if (!all(sample.ids %in% rownames(mqc))) {
      system(paste("multiqc", "-o", output.dir, output.dir))
    }
  }

  # make Table1
  bam.links <- paste0("[ bam ](./", file.path(output.dir, paste0(sample.ids, ".sorted.bam")), ")")
  bai.links <- paste0("[ bai ](./", file.path(output.dir, paste0(sample.ids, ".sorted.bam.bai")), ")")

  mqc <- read.table(file.path(output.dir, "multiqc_data", "multiqc_hisat2.txt"), sep="\t", header=TRUE)
  #mqc <- format(mqc, scientific=FALSE, big.mark=",")

  Table1 <- data.frame(SampleID=mqc$Sample,
                       Total_Reads=mqc$unpaired_total,
                       Unaligned=mqc$unpaired_aligned_none,
                       Aligned_Once=mqc$unpaired_aligned_one,
                       Aligned_Multiple=mqc$unpaired_aligned_multi,
                       Alignment_Rate=paste0(mqc$overall_alignment_rate, "%"),
                       Files=paste0(bam.links, ", ", bai.links))

  Figure1 <- function() {
    tab <- Table1[ , c(1, 3:5) ]
    #tab[,-1] <- t(t(tab[,-1]) / colSums(tab[,-1]))

    mdf <- melt(tab)

    plot(ggplot(mdf, aes(x=SampleID, weight=value, fill=variable)) +
           geom_bar() +
           ylab("Number of Reads") +
           xlab("") +
           theme(axis.text.x = element_text(angle=45, hjust=1)) +
           scale_fill_brewer(palette="Set3"))

  }

  # Setup the output chunks
  chunks <- list(
    fig1 = figure_chunk(
      title = "Barchart of mapped reads",
      description = "Barchart showing number of reads that failed to map, mapped once or mapped in multiple locations.",
      fun = Figure1,
      width= 8,
      height = 4
    ),
    tab1 = table_chunk(
      title="Summary of mapped reads",
      description="Summary of mapped reads showing the number of processed reads, number of reads that failed to align, aligned only once or aligned in multiple locations.",
      dataframe=Table1,
      collapsed = TRUE
    )
  )

  return(list(title=title,
              chunks=chunks,
              par=mget(names(formals()),sys.frame(sys.nframe())),
              description=MODULE.DESCRIPTION))
}
