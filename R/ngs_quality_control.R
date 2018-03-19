#' NGS Quality Control
#'
#' This function performs quality control of NGS datasets using FASTQC and multiqc. FastQC is run on each sample...
#'
#' @param title title of the section
#' @param sample.pathnames samples to process.
#' @param fastqc.command command to run fastqc.
#' @param output.dir output directory for results.
#' @param sample.ids optional names for the reports.
#' @param MAX.THREADS maximum number of concurrent threads to use.
#'
#' @return A list containing the following elements:
#'
#' chunks: fig1, tab1 and file1.
#'
#' @export
#'
#' @examples
ngs_quality_control <- function(title, sample.pathnames, fastqc.command, output.dir, sample.ids = NULL, MAX.THREADS = 8) {
  MODULE.DESCRIPTION <-
    "This step performs quality control of sequencing data using
  [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and
  [MultiQC](http://multiqc.info/)."

  # make a helper dataframe
  reports <- data.frame(sample.pathname = sample.pathnames,
                        sample.filename = basename(sample.pathnames))

  # set names of reports
  if (!is.null(sample.ids)) {
    reports$target.basenames <- sample.ids
  } else {
    reports$target.basenames <- gsub("\\.gz", "", reports$sample.filename)
    reports$target.basenames <- gsub("\\.fa|\\.fq|\\.fasta|\\.fastq|\\.txt", "", reports$target.basenames)
  }

  # make symbolic links to samples, allowing us to rename the samples
  tmpdir <- tempdir()
  reports$renamed.sample <- paste0(file.path(tmpdir, reports$target.basenames), ".", tools::file_ext(reports$sample.pathname))
  system(paste("ln -f -s", file.path(getwd(), reports$sample.pathname), reports$renamed.sample, collapse="\n"))

  reports$target.zip <- file.path(output.dir, paste0(reports$target.basenames, "_fastqc.zip"))
  reports$target.zip.link <- paste0("[ Zip ](", reports$target.zip, ")")
  reports$target.html <- file.path(output.dir, paste0(reports$target.basenames, "_fastqc.html"))
  reports$target.html.link <- paste0("[ HTML ](", reports$target.html, ")")

  reports$target.exists <- file.exists(file.path(reports$target.zip))

  reports$command <- paste(fastqc.command,
                           "-o", file.path(output.dir),
                           file.path(reports$renamed.sample))

  # Create directory to store FastQC files
  dir.create(file.path(output.dir), recursive = TRUE, showWarnings = FALSE)

  # Run FastQC on all samples (except if already done)
  w <- which(reports$target.exists == FALSE)
  res <- parallel::mclapply(reports$command[w], function(x) {
    system(x)
  }, mc.cores = MAX.THREADS)

  # Run MultiQC on results
  multiqc.file <- file.path(output.dir, "multiqc_report.html")

  if (!file.exists(multiqc.file)) {
    system(paste("multiqc",
                 "-o", file.path(output.dir),
                 file.path(output.dir)))
  }

  mqc <- read.table(file.path(output.dir, "multiqc_data", "multiqc_fastqc.txt"), sep="\t", header=TRUE)

  # if not all fastqc reports are in the multiqc table, redo it
  if (!all(reports$target.basenames %in% mqc$Sample)) {
    file.remove(file.path(output.dir, "multiqc_report.html"))
    unlink(file.path(output.dir, "multiqc_data"), recursive = TRUE, force = TRUE)

    system(paste("multiqc",
                 "-o", file.path(output.dir),
                 file.path(output.dir)))

    mqc <- read.table(file.path(output.dir, "multiqc_data", "multiqc_fastqc.txt"), sep="\t", header=TRUE)
  }

  # Table1
  m <-  match(reports$target.basenames, mqc$Sample)
  Table1 <- data.frame(SampleID = reports$target.basenames,
                       Filename = reports$sample.pathname,
                       Sequences = format(mqc$Total.Sequences[ m ], scientific = FALSE, big.mark = ","),
                       Sequence_Length = as.character(mqc$Sequence.length[ m ]),
                       Poor_Quality = format(mqc$Sequences.flagged.as.poor.quality[ m ], scientific = FALSE, big.mark = ","),
                       FastQC_HTML = reports$target.html.link,
                       FastQC_Zip = reports$target.zip.link)

  # Figure1
  make.fig1 <- function() {
    require(superheat)

    num <- c("pass"=1, "warn"=0, "fail"=-1)

    m <-  match(reports$target.basenames, mqc$Sample)
    vals <- mqc[ m, c(2, 4, 7, 8, 9, 12, 13, 15, 16, 17, 21, 22) ]
    vals <- apply(vals, 2, function(x) num[ match(x, names(num)) ])
    rownames(vals) <- mqc$Sample[ m ]
    colnames(vals) <- gsub("_", " ", colnames(vals))

    vals <- vals[ rev(1:nrow(vals)), ]

    superheat(X = vals,
              bottom.label.text.angle = 90,
              grid.hline.col = "white",
              grid.vline.col = "white",
              grid.hline.size = 2,
              grid.vline.size = 2,
              heat.pal = c("red4", "orange", "darkgreen"),
              left.label.size = 0.6,
              left.label.col = "white",
              left.label.text.alignment = "right",
              bottom.label.size = 1,
              bottom.label.col = "white",
              bottom.label.text.alignment = "right",
              heat.pal.values = c(0, 0.5, 1),
              legend = FALSE)
  }

  # Setup the output chunks
  chunks <- list(
    fig1 = UBIquitous::figure_chunk(
      fun = make.fig1,
      title="Summary of sample quality control",
      description='Summary figure of quality control results indicating potential biases in the samples. Legend: <span style="color:green">Pass</span>, <span style="color:orange">Warn</span>, <span style="color:red">Fail</span>.',
      width = 4 + 0.4 * nrow(reports)
    ),
    tab1 = UBIquitous::table_chunk(
      dataframe = Table1,
      title="FastQC analysis reports for all samples",
      description="Summary table of quality control results displaying basic sample statistics."
    ),
    file1 = UBIquitous::file_chunk(
      uri=file.path(output.dir, "multiqc_report.html"),
      title="MultiQC aggregated report",
      description="This report aggregates the above FastQC analyses into a single report allowing for quick comparison of dataset characteristics."
    )
  )

  par <- UBIquitous::extract_parameters()

  return(list(title=title,
              chunks=chunks,
              par=par,
              description=MODULE.DESCRIPTION))
}

