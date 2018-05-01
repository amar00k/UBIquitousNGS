# fig_clustering <- figure_chunk(
#   title="Clustering of sample to sample correlations",
#   description='Hierarchical clustering of samples based on Pearson (left) and Spearman correlations (right). The log10 of raw counts + 1 is used as input to the clustering algorithm.',
#   width = max(8, 0.5 * ncol(counts)),
#   height= max(4, 0.25 * ncol(counts)),
#   pdf.filename = file.path(output.dir, "Correlation.pdf"),
#   fun = function() {
#     log_counts <- log10(counts + 1)
#
#     par(mfrow=c(1,2))
#     h1 <- Heatmap(cor(log_counts), name="Pearson", column_title = "Pearson")
#     h2 <- Heatmap(cor(log_counts, method="spearman"), name="Spearman", column_title = "Spearman")
#
#     draw(h1 + h2)
#   }
# )


#' NGS counting with featureCounts.
#'
#' @param bam.pathnames
#' @param sample.ids
#' @param feature.counts.command
#' @param gene.annotation
#' @param output.dir
#' @param feature.counts.options
#' @param MAX.THREADS
#' @param title
#'
#' @return
#' @export
#'
#' @import ggplot2
#' @import ComplexHeatmap
#'
#' @examples
ngs_feature_counts <- function(bam.pathnames,
                               sample.ids,
                               feature.counts.command,
                               gene.annotation,
                               output.dir,
                               feature.counts.options = "-t exon -g gene_id",
                               MAX.THREADS=8,
                               title="") {

  MODULE.DESCRIPTION <- "This step quantifies gene expression by counting the number of
  mapped reads overlapping with each annotated gene using [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)."

  counts.file <- file.path(output.dir, "feature-counts.txt")

  counts.cmd <- paste(feature.counts.command,
                      "-a", gene.annotation,
                      feature.counts.options,
                      "-T", MAX.THREADS,
                      "-o", counts.file,
                      paste(bam.pathnames, collapse=" "))

  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

  log.basename <- file.path(output.dir, "ngs_feature_counts")

  # TODO: Better check if counts file include all samples.
  if (!file.exists(counts.file)) {
    run_system(counts.cmd, log.basename)
  }

  # Table: Summary of read assignment to features
  Table1 <- read.table(file.path(output.dir, "feature-counts.txt.summary"), sep="\t", header=TRUE)
  colnames(Table1)[-1] <- as.character(sample.ids)

  # simplify the feature counts table and save a new counts file
  fcounts <- read.table(counts.file, sep="\t", header=TRUE)

  counts <- fcounts[, -(1:6), drop=FALSE ]
  colnames(counts) <- as.character(sample.ids) # WARNING: is this always valid?
  rownames(counts) <- as.character(fcounts$Geneid)

  write.table(counts, file.path(output.dir, "counts.txt"), sep="\t")

  # figure 1: read assignment to features
  make.fig1 <- function() {
    w <-  which(rowSums(Table1[,-1, drop=FALSE]) > 0)

    t1 <- Table1[ w, -1, drop=FALSE ]
    status <- as.character(Table1$Status[w])

    #t1 <- t(t(t1) / colSums(t1))
    t1 <- data.frame(Status = status, t1, check.names = FALSE)
    t1$Status <- factor(t1$Status)

    mdf <- melt(t1, id.vars = "Status")
    colnames(mdf) <- c("Status", "SampleID", "Count")

    plot(ggplot(mdf, aes(x=SampleID, weight=Count, fill=Status)) +
      geom_bar() +
      ylab("Number of reads") +
      xlab("") +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      scale_fill_brewer(palette="Set3"))
  }

  # figure 2: Distribution of raw counts
  make.fig2 <- function() {
    mdf <- reshape2::melt(counts, measure.vars = colnames(counts))
    colnames(mdf) <- c("SampleID", "Counts")

    mdf <- mdf[ which(mdf$Counts > 0), ]

    p1 <- ggplot(mdf, aes(x=Counts, col=SampleID)) +
      geom_density() +
      scale_x_log10() +
      theme(legend.position = "none")

    p2 <- ggplot(mdf, aes(y=Counts, x=SampleID, col=SampleID)) +
      geom_boxplot() +
      scale_y_log10() +
      xlab("") +
      theme(axis.text.x = element_text(angle=45, hjust=1),
            legend.position = "none")

    plot(cowplot::plot_grid(p2, p1, ncol=2, labels = c("A", "B"), rel_widths = c(1,1)))
  }

  fig_clustering <- figure_chunk(
    title="Clustering of sample to sample correlations",
    description='Hierarchical clustering of samples based on Pearson (left) and Spearman correlations (right). The log10 of raw counts + 1 is used as input to the clustering algorithm.',
    width = max(8, 0.5 * ncol(counts)),
    height= max(4, 0.25 * ncol(counts)),
    pdf.filename = file.path(output.dir, "Correlation.pdf"),
    fun = function() {
      log_counts <- log10(counts + 1)

      par(mfrow=c(1,2))
      h1 <- Heatmap(cor(log_counts), name="Pearson", column_title = "Pearson")
      h2 <- Heatmap(cor(log_counts, method="spearman"), name="Spearman", column_title = "Spearman")

      draw(h1 + h2)
    }
  )

  # Setup the output chunks
  chunks <- list(
    commands = commands_chunk(
      commands.filename = paste0(log.basename, ".commands.sh"),
      stdout.filename = paste0(log.basename, ".stdout.txt"),
      stderr.filename = paste0(log.basename, ".stderr.txt"),
      title = "Feature counts commands",
      description="Commands used to count reads in features.",
      collapsed=TRUE),
    fig1 = figure_chunk(
      fun = make.fig1,
      title="Barchart of read assignment to features",
      description='This figure summarizes read assignment to genomic features.',
      width = 8,
      height= 4,
      pdf.filename = file.path(output.dir, "Features.pdf")
    ),
    tab1 = table_chunk(
      title="Summary of read assignment to features",
      description="This table summarizes read assignment to genomic features.",
      scrollX = TRUE,
      dataframe=Table1,
      collapsed = TRUE
    ),
    fig2 = figure_chunk(
      fun = make.fig2,
      title="Distributions of raw counts",
      description='Distributions of log10 raw counts with zeros removed. **(A)** Boxplots; **(B)** Density plots.',
      width = 8,
      height= 4,
      pdf.filename = file.path(output.dir, "Distribution.pdf")
    ),
    fig3 = fig_clustering
  )

  return(list(title=title,
              chunks=chunks,
              par=mget(names(formals()), sys.frame(sys.nframe())),
              description=MODULE.DESCRIPTION,
              counts=counts))
}
