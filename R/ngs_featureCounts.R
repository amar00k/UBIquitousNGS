library(UBIquitous)

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
#' @examples
ngs_feature_counts <- function(bam.pathnames,
                               sample.ids,
                               feature.counts.command,
                               gene.annotation,
                               output.dir,
                               feature.counts.options = "-t exon -g gene_id",
                               MAX.THREADS=8,
                               title="") {
  require(ggplot2)
  require(reshape2)
  require(cowplot)

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

  # TODO: Better check if counts file include all samples.
  if (!file.exists(counts.file)) {
    system(counts.cmd)
  }

  # Table: Summary of read assignment to features
  Table1 <- read.table(file.path(output.dir, "feature-counts.txt.summary"), sep="\t", header=TRUE)
  colnames(Table1)[-1] <- as.character(sample.ids)

  # simplify the feature counts table and save a new counts file
  fcounts <- read.table(counts.file, sep="\t", header=TRUE)

  counts <- fcounts[, -(1:6) ]
  colnames(counts) <- as.character(sample.ids) # WARNING: is this always valid?
  rownames(counts) <- as.character(fcounts$Geneid)

  write.table(counts, file.path(output.dir, "counts.txt"), sep="\t")

  # figure 1: read assignment to features
  make.fig1 <- function() {
    t1 <- Table1[ which(rowSums(Table1[,-1]) > 0), ]
    t1$Status <- droplevels(t1$Status)
    t1[, -1] <- t(t(t1[, -1]) / colSums(t1[, -1]))

    mdf <- melt(t1)
    colnames(mdf) <- c("Status", "SampleID", "Count")

    plot(ggplot(mdf, aes(x=SampleID, weight=Count, fill=Status)) +
      geom_bar() +
      ylab("Relative Abundance") +
      xlab("") +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      scale_fill_brewer(palette="Set3"))
  }

  # figure 2: Distribution of raw counts
  make.fig2 <- function() {
    mdf <- melt(counts)
    colnames(mdf) <- c("SampleID", "Counts")

    p1 <- ggplot(mdf, aes(x=Counts+1, col=SampleID)) +
      geom_density() +
      xlab("Total Counts") +
      scale_x_log10() +
      theme(legend.position = "none")

    p2 <- ggplot(mdf, aes(y=Counts+1, x=SampleID, col=SampleID)) +
      geom_boxplot() +
      scale_y_log10() +
      xlab("") +
      theme(axis.text.x = element_text(angle=45, hjust=1),
            legend.position = "none")

    plot(cowplot::plot_grid(p2, p1, ncol=2, labels = c("A", "B"), rel_widths = c(1,1)))
  }

  # figure 3: Clustering of samples
  make.fig3 <- function() {
    require(ComplexHeatmap)

    par(mfrow=c(1,2))
    h1 <- Heatmap(cor(counts), name="Pearson", column_title = "Pearson")
    h2 <- Heatmap(cor(counts, method="spearman"), name="Spearman", column_title = "Spearman")

    draw(h1 + h2)
  }

  # Setup the output chunks
  chunks <- list(
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
      description='Distributions of raw counts. **(A)** Boxplots; **(B)** Density plots.',
      width = 8,
      height= 4,
      pdf.filename = file.path(output.dir, "Distribution.pdf")
    ),
    fig3 = figure_chunk(
      fun = make.fig3,
      title="Clustering of sample to sample correlations",
      description='Hierarchical clustering of samples based on Pearson (left) and Spearman correlations (right).',
      width = 8,
      height= 4,
      pdf.filename = file.path(output.dir, "Correlation.pdf")
    )
  )

  return(list(title=title,
              chunks=chunks,
              par=mget(names(formals()), sys.frame(sys.nframe())),
              description=MODULE.DESCRIPTION))
}
