
#' NGS multiqc
#'
#' @param input.dir
#' @param output.dir
#' @param multiqc.command
#' @param log.basename
#' @param overwrite
#'
#' @return
#' @export
#'
#' @examples
ngs_multiqc <- function(input.dir,
                        output.dir,
                        multiqc.command = "multiqc",
                        multiqc.options = "",
                        log.basename = NULL,
                        overwrite = FALSE) {
  multiqc.file <- file.path(output.dir, "multiqc_report.html")

  if (is.null(log.basename)) {
    log.basename <- file.path(output.dir, "multiqc_system")
  }

  # do multiqc if there have been modifications in the input folder
  do.command = TRUE

  if (file.exists(multiqc.file)) {
    last.change <- max(file.info(list.files(input.dir, full.names = TRUE))$mtime)

    if (last.change <= file.info(multiqc.file)$mtime) {
      do.command = FALSE
    }
  }

  command <- paste(multiqc.command,
                   multiqc.options,
                   "-f",
                   "-o", file.path(output.dir),
                   file.path(input.dir))

  if (do.command == TRUE) {
    run_system(command, log.basename = log.basename)
  }

  # check for reports, and make chunks
  chunks <- list()
  # script_command = script_chunk(
  #   commands = command,
  #   title = "MultiQC commands",
  #   description = "Commands used to generate a MultiQC report.",
  #   collapsed = TRUE
  # ))

  fastqc.report <- file.path(output.dir, "multiqc_data", "multiqc_fastqc.txt")
  trimmomatic.report <- file.path(output.dir, "multiqc_data", "multiqc_trimmomatic.txt")
  hisat2.report <- file.path(output.dir, "multiqc_data", "multiqc_hisat2.txt")

  # Handle fastqc report
  if (file.exists(fastqc.report)) {
    mqc <- read.table(fastqc.report, sep="\t", header=TRUE, stringsAsFactors = FALSE)

    chunks$figure_fastqc <- figure_chunk(
      title="Summary of sample quality control",
      description='Summary figure of quality control results indicating potential biases or problems in the samples. Legend: <span style="color:green">Pass</span>, <span style="color:orange">Warn</span>, <span style="color:red">Fail</span>.',
      width = 6,
      height = 2 + 0.15 * nrow(mqc),
      fun = function() {
        vals <- mqc[ , c(2, 4, 7, 8, 9, 12, 13, 15, 16, 17, 21, 22) ]
        rownames(vals) <- mqc$Sample

        mdf <- reshape2::melt(as.matrix(vals),
                              measure.vars = colnames(vals),
                              varnames = c("SampleID", "Test"))
        mdf$SampleID <- factor(mdf$SampleID, levels=rev(levels(mdf$SampleID)))

        ggplot(mdf, aes(x=Test, y=SampleID)) +
          geom_tile(aes(fill=value), col="white", lwd=1.5) +
          scale_fill_manual(values = c(fail="red4", pass="darkgreen", warn = "orange")) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      })

    df <- data.frame(SampleID = mqc$Sample,
                     Sequences = format(mqc$Total.Sequences, scientific = FALSE, big.mark = ","),
                     Lengths = as.character(mqc$Sequence.length),
                     Avg_Length = format(mqc$avg_sequence_length, digits = 1, nsmall = 1),
                     Poor_Quality = format(mqc$Sequences.flagged.as.poor.quality, scientific = FALSE, big.mark = ","),
                     HTML = paste0("[HTML](", file.path(input.dir, paste0(mqc$Sample, "_fastqc.html")), ")"),
                     Zip = paste0("[Zip](", file.path(input.dir, paste0(mqc$Sample, "_fastqc.zip")), ")"))

    chunks$table_fastqc <- table_chunk(
      dataframe = df,
      title="FastQC analysis reports for all samples",
      description="Summary table of quality control results displaying basic sample statistics.",
      collapsed = TRUE,
      scrollX = FALSE)
  }

  # Handle Trimmomatic report
  # TODO: Handle single-end data!
  if (file.exists(trimmomatic.report)) {
    mqc <- read.table(trimmomatic.report, sep="\t", header=TRUE, stringsAsFactors = FALSE)

    chunks$figure_trimmomatic <- figure_chunk(
      title="Summary of Trimmomatic results",
      description='Summary of trimming procedure showing the number of surviving and dropped reads.',
      width = 4 + nrow(mqc)*0.25,
      height = 5,
      fun = function() {
        df <- mqc[, c("Sample", "surviving", "forward_only_surviving", "reverse_only_surviving", "dropped") ]
        mdf <- reshape2::melt(df, id.vars = "Sample")

        ggplot(mdf, aes(x=Sample, weight=value, fill=variable)) +
          geom_bar() +
          theme(axis.text.x = element_text(angle=45, hjust=1)) +
          scale_fill_brewer(palette="Set3") +
          ylab("Number of reads")
      })
  }

  # Handle hisat2 report
  if (file.exists(hisat2.report)) {
    mqc <- read.table(hisat2.report, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    colnames(mqc)[1] <- "SampleID"

    chunks$figure_hisat2 <- figure_chunk(
      title = "Barchart of mapped reads",
      description = "Barchart showing number of reads that failed to map, mapped once or mapped in multiple locations.",
      width= 8,
      height = 4,
      fun = function() {
        tab <- mqc
        tab <- tab[, which(colnames(tab) %in% c("SampleID", "paired_aligned_none", "paired_aligned_discord_one", "unpaired_aligned_one",	"unpaired_aligned_none", "unpaired_aligned_multi", "paired_aligned_multi", "paired_aligned_one")) ]
        tab$SampleID <- factor(tab$SampleID, levels=tab$SampleID)
        #tab[,-1] <- t(t(tab[,-1]) / colSums(tab[,-1]))

        mdf <- reshape2::melt(tab, id.vars="SampleID")

        plot(ggplot(mdf, aes(x=SampleID, weight=value, fill=variable)) +
               geom_bar() +
               ylab("Number of Reads") +
               xlab("") +
               theme(axis.text.x = element_text(angle=45, hjust=1)) +
               scale_fill_brewer(palette="Set3"))

      }
    )

    chunks$table_hisat2 <- table_chunk(
      title="Summary of mapped reads",
      description="Summary of mapped reads showing the number of processed reads, number of reads that failed to align, aligned only once or aligned in multiple locations.",
      dataframe=mqc,
      collapsed = TRUE,
      scrollX = TRUE
    )
  }

  return(list(chunks = chunks))
}

#' NGS FastQC
#'
#' This function performs quality control of NGS datasets using FASTQC and multiqc. FastQC is run on each sample...
#'
#' @import UBIquitous
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
ngs_fastqc <- function(sample.pathnames,
                       output.dir,
                       sample.ids = NULL,
                       fastqc.command = "fastqc",
                       multiqc.command = "multiqc",
                       MAX.THREADS = 8,
                       title="") {
  MODULE.DESCRIPTION <-
    "We performed quality control of the sequencing data using
  [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and
  [MultiQC](http://multiqc.info/)."

  # make a helper dataframe
  reports <- data.frame(sample.pathname = sample.pathnames,
                        sample.filename = basename(sample.pathnames))

  reports$sample.extension <- ifelse(tools::file_ext(reports$sample.pathname) == "gz",
                                     paste0(tools::file_ext(gsub("\\.gz$", "", reports$sample.pathname)), ".gz"),
                                     tools::file_ext(reports$sample.pathname))

  # set names of reports
  if (!is.null(sample.ids)) {
    reports$target.basenames <- sample.ids
  } else {
    reports$target.basenames <- gsub("\\.gz$", "", reports$sample.filename)
    reports$target.basenames <- gsub("\\.fa$|\\.fq$|\\.fasta$|\\.fastq$|\\.txt$", "", reports$target.basenames)
  }

  # cleanup files that don't exist (shoud emit a waring somehow)
  reports <- reports[ which(file.exists(file.path(reports$sample.pathname))), ]

  # make symbolic links to samples, allowing us to rename the samples
  sample.dir <- file.path(output.dir, "samples_fastqc")
  reports$renamed.sample <- paste0(file.path(sample.dir, reports$target.basenames), ".", reports$sample.extension)

  reports$target.zip <- file.path(output.dir, paste0(reports$target.basenames, "_fastqc.zip"))
  reports$target.html <- file.path(output.dir, paste0(reports$target.basenames, "_fastqc.html"))

  reports$rename.command <- paste("ln -f", reports$sample.pathname, reports$renamed.sample)
  reports$fastqc.command <- paste(fastqc.command, "-o", output.dir, reports$renamed.sample)

  reports$target.exists <- file.exists(file.path(reports$target.zip))

  # run the commands
  log.basename <- file.path(output.dir, "fastqc_system")

  if (!dir.exists(output.dir)) {
    run_system(paste("mkdir -p", output.dir), log.basename)
  }

  if (!dir.exists(sample.dir)) {
    run_system(paste("mkdir -p", sample.dir), log.basename)
  }

  w <- which(reports$target.exists == FALSE)

  # TODO: check why sometimes this gets stuck
  # It's because FastQC hangs with reads that are too small or zero length
  commands <- reports$rename.command[w]
  run_system(commands, log.basename)

  commands <- reports$fastqc.command[w]
  run_system_parallel(commands, log.basename, max.cores = MAX.THREADS)

  # if (length(w) > 0) {
  #   commands <- reports$rename.command[w]
  #   run_system(commands, log.basename)
  #
  #   commands <- paste(fastqc.command,
  #                     "-o", output.dir,
  #                     "-t", MAX.THREADS,
  #                     paste(reports$renamed.sample[w], collapse=" "))
  #
  #   run_system(commands, log.basename)
  # }

  # Run MultiQC on results
  multiqc <- ngs_multiqc(input.dir = file.path(output.dir),
                         output.dir = file.path(output.dir, "multiqc"),
                         multiqc.command = multiqc.command,
                         log.basename = log.basename)

  # Setup the output chunks
  chunks <- list(
    commands_fastqc = commands_chunk(
      commands.filename = paste0(log.basename, ".commands.sh"),
      stdout.filename = paste0(log.basename, ".stdout.txt"),
      stderr.filename = paste0(log.basename, ".stderr.txt"),
      title = "Quality control commands",
      description="Commands executed to generate FastQC and MultiQC reports for all samples.",
      collapsed=TRUE),
    figure_fastqc = multiqc$chunks$figure_fastqc,
    table_fastqc = multiqc$chunks$table_fastqc
    # file1 = file_chunk(
    #   uri=file.path(output.dir, "multiqc_report.html"),
    #   title="MultiQC aggregated report",
    #   description="This report aggregates the above FastQC analyses into a single report allowing for quick comparison of dataset characteristics."
    # )
  )

  par <- extract_parameters()

  return(list(title=title,
              chunks=chunks,
              par=par,
              description=MODULE.DESCRIPTION))
}


















#'
#'
#' #' NGS Quality Control
#' #'
#' #' This function performs quality control of NGS datasets using FASTQC and multiqc. FastQC is run on each sample...
#' #'
#' #' @import UBIquitous
#' #'
#' #' @param title title of the section
#' #' @param sample.pathnames samples to process.
#' #' @param fastqc.command command to run fastqc.
#' #' @param output.dir output directory for results.
#' #' @param sample.ids optional names for the reports.
#' #' @param MAX.THREADS maximum number of concurrent threads to use.
#' #'
#' #' @return A list containing the following elements:
#' #'
#' #' chunks: fig1, tab1 and file1.
#' #'
#' #' @export
#' #'
#' #' @examples
#' ngs_quality_control <- function(sample.pathnames, fastqc.command, output.dir, sample.ids = NULL, MAX.THREADS = 8, title="") {
#'   MODULE.DESCRIPTION <-
#'     "This step performs quality control of sequencing data using
#'   [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and
#'   [MultiQC](http://multiqc.info/)."
#'
#'   # make a helper dataframe
#'   reports <- data.frame(sample.pathname = sample.pathnames,
#'                         sample.filename = basename(sample.pathnames))
#'
#'   # set names of reports
#'   if (!is.null(sample.ids)) {
#'     reports$target.basenames <- sample.ids
#'   } else {
#'     reports$target.basenames <- gsub("\\.gz", "", reports$sample.filename)
#'     reports$target.basenames <- gsub("\\.fa|\\.fq|\\.fasta|\\.fastq|\\.txt", "", reports$target.basenames)
#'   }
#'
#'   # make symbolic links to samples, allowing us to rename the samples
#'   tmpdir <- tempdir()
#'   reports$renamed.sample <- paste0(file.path(tmpdir, reports$target.basenames), ".", tools::file_ext(reports$sample.pathname))
#'   system(paste("ln -f -s", file.path(getwd(), reports$sample.pathname), reports$renamed.sample, collapse="\n"))
#'
#'   reports$target.zip <- file.path(output.dir, paste0(reports$target.basenames, "_fastqc.zip"))
#'   reports$target.zip.link <- paste0("[ Zip ](", reports$target.zip, ")")
#'   reports$target.html <- file.path(output.dir, paste0(reports$target.basenames, "_fastqc.html"))
#'   reports$target.html.link <- paste0("[ HTML ](", reports$target.html, ")")
#'
#'   reports$target.exists <- file.exists(file.path(reports$target.zip))
#'
#'   reports$command <- paste(fastqc.command,
#'                            "-o", file.path(output.dir),
#'                            file.path(reports$renamed.sample))
#'
#'   # Create directory to store FastQC files
#'   dir.create(file.path(output.dir), recursive = TRUE, showWarnings = FALSE)
#'
#'   # Run FastQC on all samples (except if already done)
#'   w <- which(reports$target.exists == FALSE)
#'   res <- parallel::mclapply(reports$command[w], function(x) {
#'     system(x)
#'   }, mc.cores = MAX.THREADS)
#'
#'   # Run MultiQC on results
#'   multiqc.file <- file.path(output.dir, "multiqc_report.html")
#'
#'   if (!file.exists(multiqc.file)) {
#'     system(paste("multiqc",
#'                  "-o", file.path(output.dir),
#'                  file.path(output.dir)))
#'   }
#'
#'   mqc <- read.table(file.path(output.dir, "multiqc_data", "multiqc_fastqc.txt"), sep="\t", header=TRUE)
#'
#'   # if not all fastqc reports are in the multiqc table, redo it
#'   if (!all(reports$target.basenames %in% mqc$Sample)) {
#'     file.remove(file.path(output.dir, "multiqc_report.html"))
#'     unlink(file.path(output.dir, "multiqc_data"), recursive = TRUE, force = TRUE)
#'
#'     system(paste("multiqc",
#'                  "-o", file.path(output.dir),
#'                  file.path(output.dir)))
#'
#'     mqc <- read.table(file.path(output.dir, "multiqc_data", "multiqc_fastqc.txt"), sep="\t", header=TRUE)
#'   }
#'
#'   # Table1
#'   m <-  match(reports$target.basenames, mqc$Sample)
#'   Table1 <- data.frame(SampleID = reports$target.basenames,
#'                        Filename = reports$sample.pathname,
#'                        Sequences = format(mqc$Total.Sequences[ m ], scientific = FALSE, big.mark = ","),
#'                        Sequence_Length = as.character(mqc$Sequence.length[ m ]),
#'                        Poor_Quality = format(mqc$Sequences.flagged.as.poor.quality[ m ], scientific = FALSE, big.mark = ","),
#'                        FastQC_HTML = reports$target.html.link,
#'                        FastQC_Zip = reports$target.zip.link)
#'
#'   # Figure1
#'   make.fig1 <- function() {
#'     num <- c("pass"=1, "warn"=0, "fail"=-1)
#'
#'     m <-  match(reports$target.basenames, mqc$Sample)
#'     vals <- mqc[ m, c(2, 4, 7, 8, 9, 12, 13, 15, 16, 17, 21, 22) ]
#'     vals <- apply(vals, 2, function(x) num[ match(x, names(num)) ])
#'     rownames(vals) <- mqc$Sample[ m ]
#'     colnames(vals) <- gsub("_", " ", colnames(vals))
#'
#'     vals <- vals[ rev(1:nrow(vals)), ]
#'
#'     superheat::superheat(X = vals,
#'               bottom.label.text.angle = 90,
#'               grid.hline.col = "white",
#'               grid.vline.col = "white",
#'               grid.hline.size = 2,
#'               grid.vline.size = 2,
#'               heat.pal = c("red4", "orange", "darkgreen"),
#'               left.label.size = 0.6,
#'               left.label.col = "white",
#'               left.label.text.alignment = "right",
#'               bottom.label.size = 0.5,
#'               bottom.label.col = "white",
#'               bottom.label.text.alignment = "right",
#'               heat.pal.values = c(0, 0.5, 1),
#'               legend = FALSE)
#'   }
#'
#'   # Setup the output chunks
#'   chunks <- list(
#'     fig1 = figure_chunk(
#'       fun = make.fig1,
#'       title="Summary of sample quality control",
#'       description='Summary figure of quality control results indicating potential biases in the samples. Legend: <span style="color:green">Pass</span>, <span style="color:orange">Warn</span>, <span style="color:red">Fail</span>.',
#'       width = 6,
#'       height = 4 + 0.2 * nrow(reports)
#'     ),
#'     tab1 = table_chunk(
#'       dataframe = Table1,
#'       title="FastQC analysis reports for all samples",
#'       description="Summary table of quality control results displaying basic sample statistics."
#'     ),
#'     file1 = file_chunk(
#'       uri=file.path(output.dir, "multiqc_report.html"),
#'       title="MultiQC aggregated report",
#'       description="This report aggregates the above FastQC analyses into a single report allowing for quick comparison of dataset characteristics."
#'     )
#'   )
#'
#'   par <- extract_parameters()
#'
#'   return(list(title=title,
#'               chunks=chunks,
#'               par=par,
#'               description=MODULE.DESCRIPTION))
#' }
#'
