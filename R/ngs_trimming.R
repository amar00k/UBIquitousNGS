#' NGS Trimming
#'
#' Performs trimming of fastq files.
#'
#' @import UBIquitous
#'
#' @param sample.pathnames samples to process.
#' @param cutadapt.options options to pass to cutadapt command line.
#' @param cutadapt.command command to run cutadapt.
#' @param output.dir directory to store processed samples.
#' @param sample.ids basenames to use for output files.
#' @param MAX.THREADS maximum number of concurrent threads to run.
#'
#' @return
#' @export
#'
#' @examples
ngs_trimming <- function(sample.pathnames,
                         cutadapt.options,
                         cutadapt.command,
                         output.dir,
                         sample.ids = NULL,
                         MAX.THREADS=8,
                         title="") {

  MODULE.DESCRIPTION <- "This step performs trimming of the raw reads using cutadapt
  in order to remove adapter sequences still present in the reads as well as removing
  low-quality bases from the ends of the reads."

  data <- data.frame(sample.pathname = sample.pathnames,
                     sample.filename = basename(sample.pathnames),
                     stringsAsFactors = FALSE)

  if (!is.null(sample.ids)) {
    target.extension <- ifelse(tools::file_ext(data$sample.filename) == "gz",
                               paste0(tools::file_ext(gsub("\\.gz", "", data$sample.filename)), ".gz"),
                               tools::file_ext(data$sample.filename))

    data$target.filename <- file.path(output.dir, paste0(sample.ids, ".", target.extension))
  } else {
    data$target.filename <- file.path(output.dir, basename(data$sample.pathname))
  }

  data$target.basenames <- gsub("\\.gz", "", basename(data$target.filename))
  data$target.basenames <- gsub("\\.fa|\\.fq|\\.fasta|\\.fastq|\\.txt", "", data$target.basenames)

  data$target.log <- paste0(data$target.filename, ".log")
  data$target.log.link <- paste0("[ Log ](", data$target.log, ")")

  data$target.exists <- file.exists(file.path(data$target.filename))

  data$command <- paste(cutadapt.command,
                        cutadapt.options,
                        data$sample.pathname,
                        "-o", data$target.filename,
                        ">", data$target.log)

  # create dir and do trimming
  dir.create(file.path(output.dir), recursive = TRUE, showWarnings = FALSE)

  w <- which(data$target.exists == FALSE)
  res <- parallel::mclapply(data$command[w], function(x) {
    system(x)
  }, mc.cores = MAX.THREADS)

  # do multiqc
  multiqc.file <- file.path(output.dir, "multiqc_report.html")

  if (!file.exists(multiqc.file)) {
    system(paste("multiqc",
                 "-o", file.path(output.dir),
                 file.path(output.dir)))
  }

  mqc <- read.table(file.path(output.dir, "multiqc_data", "multiqc_cutadapt.txt"), sep="\t", header=TRUE)

  # make Table1
  m <-  match(data$target.basenames, mqc$Sample)
  Table1 <- data.frame(SampleID = data$target.basenames,
                    Filename = data$sample.pathname,
                    Processed = format(mqc$r_processed[ m ], scientific = FALSE, big.mark = ","),
                    Adapters_Trimmed = format(mqc$r_with_adapters[ m ], scientific = FALSE, big.mark = ","),
                    Quality_Trimmed = format(mqc$quality_trimmed[ m ], scientific = FALSE, big.mark = ","),
                    Percent_Trimmed = paste0(format(mqc$percent_trimmed[ m ], digits=2), "%"),
                    Log = data$target.log.link)

  # Setup the output chunks
  chunks <- list(
    tab1 = table_chunk(
      title="Summary of trimming procedure",
      description="Summary of trimming procedure showing number of processed reads, number of reads from which adapter sequences were removed and number of reads that were quality trimmed.",
      dataframe=Table1
    )
  )

  par <- UBIquitous::extract_parameters()

  return(list(title=title, chunks=chunks, par=par, description=MODULE.DESCRIPTION))
}
