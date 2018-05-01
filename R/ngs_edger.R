
#' NGS differential expression with edgeR
#'
#' @param counts.file
#' @param sample.ids
#' @param condition.factor
#' @param output.basename
#' @param output.dir
#' @param cutoff.fdr
#' @param cutoff.lfc
#'
#' @return
#' @export
#'
#' @examples
ngs_edger <- function(counts.file,
                      sample.ids,
                      condition.factor,
                      output.basename,
                      output.dir,
                      cutoff.fdr = 0.05,
                      cutoff.lfc = 1,
                      title = "") {
  require(edgeR)
  require(ggplot2)
  require(reshape2)
  require(cowplot)

  MODULE.DESCRIPTION = "This step tests for pairwise differential expression using edgeR."

  ## define edgeR functions
  plot.volcano <- function(result) {
    colors <- ifelse(result$FDR < 0.05 & abs(result$logFC) >= 1, "red", "black")
    plot(result$logFC, -log10(result$FDR), xlab="log2 FC", ylab="-log10 FDR", col=colors, pch=20, cex=0.75)
    abline(v=c(-1, 1), lty="dashed", col="grey")
    abline(h=-log10(0.05), lty="dashed", col="grey")
  }

  plot.ma <- function(result) {
    colors <- ifelse(result$FDR < 0.05 & abs(result$logFC) >= 1, "red", "black")
    plot(result$logCPM, result$logFC, xlab="Average log CPM", ylab="log2 FC", col=colors, pch=20, cex=0.75)
    abline(h=c(-1, 1), lty="dashed", col="grey")
  }

  do.edger <- function(counts, columns, condition.factor, basename, target.dir) {
    dir.create(target.dir, showWarnings = FALSE)

    output.filename <- file.path(target.dir, paste0(basename, ".differential_expression_edgeR.xls"))

    if (file.exists(output.filename))
      return(read.table(output.filename, sep="\t", header=TRUE, check.names = FALSE))

    counts <- counts[, match(columns, colnames(counts)) ]

    # remove zeros
    w <- which(rowSums(counts) > 0)
    counts <- counts[w, ]

    y <- DGEList(counts=counts, group=condition.factor)
    y <- calcNormFactors(y)

    y <- estimateDisp(y)
    et <- exactTest(y)

    topgenes <- topTags(et, n=dim(counts)[[1]])

    norm <- cpm(y)
    colnames(norm) <- paste0(colnames(norm), "_", "CPM")
    colnames(counts) <- paste0(colnames(counts), "_", "Raw")

    cts <- merge(counts, norm, by="row.names")
    colnames(cts)[1] <- "GeneID"
    result <- merge(cts, as.data.frame(topgenes), by.x="GeneID", by.y="row.names")

    result$Result <- ifelse(result$FDR < 0.05 & abs(result$logFC) >= 1, ifelse(result$logFC < 0, "Down", "Up"), "No Change")

    result <- result[ order(result$FDR), ]

    write.table(result, output.filename, sep="\t", quote=FALSE, row.names=FALSE)

    return(result)
  }

  # read the raw counts
  counts <- read.table(counts.file, sep="\t", header=TRUE, check.names=FALSE)

  # run the test
  res <- do.edger(counts, sample.ids, condition.factor, output.basename, output.dir)

  # table of differential expressed genes
  res.diff <- res[ which(res$Result != "No Change"), ]

  # Figure: MDS plot
  make.fig1 <- function() {
    colors <- c("orange", "blue")[ as.numeric(condition.factor) ]

    cts <- res[, grep("_Raw$", colnames(res)) ]
    colnames(cts) <- gsub("_Raw$", "", colnames(cts))
    y <- DGEList(counts=cts, group=condition.factor)
    y <- calcNormFactors(y)
    plotMDS(y, col=colors, cex=0.75)
  }

  # Figure: Volcano and MA plots
  make.fig2 <- function() {
    par(mfrow=c(1,2))
    plot.volcano(res)
    plot.ma(res)
  }

  # Setup the output chunks
  chunks <- list(
    fig1 = figure_chunk(
      fun = make.fig1,
      title="MDS plot",
      description=paste0('Multi-dimensional scaling (MDS) plot showing similarity between samples. Samples are colored according to condition: <span style="color:orange"> **',
                         levels(condition.factor)[1], '** in orange</span> and <span style="color:blue"> **',
                         levels(condition.factor)[2], '** in blue</span>.'),
      width = 5,
      height= 5,
      pdf.filename = file.path(output.dir, paste0(output.basename, "-MDS-plot.pdf"))
    ),
    fig2 = figure_chunk(
      fun = make.fig2,
      title="Volcano and MA plots",
      description=paste0('**(A)** Volcano plot showing the relationship between fold-change and evidence of differential expression (-log FDR). **(B)** MA-plot showing relationship between mean expression and fold-change. Differentially expressed genes are indicated in <span style="color:red"> red </span> (FDR < ', cutoff.fdr, ' and abs(log2 FC) > ', cutoff.lfc, ').'),
      width = 10,
      height= 5,
      pdf.filename = file.path(output.dir, paste0(output.basename, "-Volcano-MA-plot.pdf"))
    ),
    tab1 = table_chunk(
      dataframe = res.diff,
      scrollX = TRUE,
      title = "Differentially expressed genes",
      collapsed = TRUE
    )
  )

  return(list(title=title,
              chunks=chunks,
              par=mget(names(formals())),
              description=MODULE.DESCRIPTION,
              result.filename=file.path(output.dir, paste0(output.basename, ".differential_expression_edgeR.xls")),
              result=res))
}


