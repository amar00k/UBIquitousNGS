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
#' @import ggplot2
#'
#' @examples
ngs_hisat2 <- function(reads1.pathnames,
                       sample.ids,
                       hisat2.command,
                       hisat2.build.command,
                       genome.fasta,
                       genome.index.name,
                       output.dir,
                       multiqc.command = "multiqc",
                       MAX.THREADS=8,
                       reads2.pathnames = rep(NA, length(reads1.pathnames)),
                       title="") {
  MODULE.DESCRIPTION <- "This step maps the sequenced reads to a reference genome
  using [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)."

  # log
  log.basename <- file.path(output.dir, "ngs_hisat2")

  # build the index
  cmd.build <- paste(hisat2.build.command,
                     "-p", MAX.THREADS,
                     genome.fasta,
                     genome.index.name)

  file.check <- paste0(genome.index.name, ".1.ht2")

  # check for index and build
  res <- integer(0)
  if (!file.exists(file.check)) {
    #message <- system(paste(cmd.build, "2>&1"), intern=TRUE)
    #res <- attr(message, "status")
    run_system(cmd.build, log.basename = log.basename)
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
                     "-p", MAX.THREADS,
                     "--new-summary --summary-file", file.path(target.dir, paste0(sample.id, ".log")),
                     "-x", genome.index.name,
                     flag.reads,
                     ">", file.path(target.dir, paste0(sample.id, ".sam")))

    #paste(c(map.cmd, sam.to.bam.commands(target.dir, sample.id)), collapse="\n")
    #paste(c(map.cmd, sam.to.bam.commands(target.dir, sample.id)), collapse="\n")
  }

  target.sams <- file.path(output.dir, paste0(sample.ids, ".sam"))
  target.files <- file.path(output.dir, paste0(sample.ids, ".sorted.bam"))

  map.cmds <- unlist(lapply(1:length(sample.ids),
                            function(i) make.hisat2.commands(reads1.pathnames[i],
                                                             reads2.pathnames[i],
                                                             sample.ids[i],
                                                             output.dir)))
  convert.commands <- sapply(1:length(sample.ids),
                             function(i) paste(sam.to.bam.commands(output.dir, sample.ids[i]), collapse="\n"))

  # dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output.dir)) {
    run_system(paste("mkdir -p", output.dir), log.basename)
  }

  w <- which(!file.exists(target.files))
  wsam <- which(!file.exists(target.sams))
  wsam <- intersect(wsam, w)

  run_system(map.cmds[wsam], log.basename = log.basename)
  run_system_parallel(convert.commands[w], log.basename = log.basename, max.cores = MAX.THREADS)

  # Run MultiQC on results
  multiqc <- ngs_multiqc(input.dir = file.path(output.dir),
                         output.dir = file.path(output.dir, "multiqc"),
                         multiqc.command = multiqc.command,
                         multiqc.options = paste("-x", paste0(basename(log.basename), ".stderr.txt")), # exclude
                         log.basename = log.basename)

  # make Table1
  bam.links <- paste0("[ bam ](./", file.path(output.dir, paste0(sample.ids, ".sorted.bam")), ")")
  bai.links <- paste0("[ bai ](./", file.path(output.dir, paste0(sample.ids, ".sorted.bam.bai")), ")")

  # mqc <- read.table(file.path(output.dir, "multiqc", "multiqc_data", "multiqc_hisat2.txt"), sep="\t", header=TRUE)
  # #mqc <- format(mqc, scientific=FALSE, big.mark=",")
  #
  # Table1 <- mqc
  # colnames(Table1)[1] <- "SampleID"
  # Table1$Files <- paste0(bam.links, ", ", bai.links)
  # Table1 <- Table1[ match(sample.ids, Table1$SampleID), ]

  commands <- commands_chunk(
    commands.filename = paste0(log.basename, ".commands.sh"),
    stdout.filename = paste0(log.basename, ".stdout.txt"),
    stderr.filename = paste0(log.basename, ".stderr.txt"),
    title = "Mapping with Hisat2",
    description="Commands executed to map reads with Hisat2.",
    collapsed=TRUE)

  return(list(title=title,
              chunks=list(
                commands = commands,
                figure_hisat2 = multiqc$chunks$figure_hisat2,
                table_hisat2 = multiqc$chunks$table_hisat2
              ),
              par=mget(names(formals()),sys.frame(sys.nframe())),
              description=MODULE.DESCRIPTION,
              bam.pathnames=target.files,
              commands=map.cmds))
}
