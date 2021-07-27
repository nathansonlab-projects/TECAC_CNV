
#  pluta  7/8/21

# map of gene and CNVR overlap

# example annotation data
# annot.dat <- read.table("call_gene/interval_annotation.txt", header = TRUE, sep = "\t")

library(ggplot2)

# ----------------------------- plotGenesInCNVR ---------------------------------- #
plotGenesInCNVR <- function(cnvr.id, cnvr.dat, annot.dat)
#  input: cnvr.id (string), id of the CNVR, e.g "CNVR_190"
#         cnvr.dat (data.frame), df of CNVR information (chr and start/end); the output of plotGeneInCNVR
#         annot.dat (data.frame), from the interval_annotation file created by HandyCNV; the output of call_gene
#
#  output: p1  (ggplot object), a plot of the CNVR and overlapping genes

{
  
  cnvr.dat <- data.frame( chr = cnvr.dat$Chr[ cnvr.dat$CNVR_ID == cnvr.id ],
                          start = cnvr.dat$Start[ cnvr.dat$CNVR_ID == cnvr.id ],
                          end = cnvr.dat$End[ cnvr.dat$CNVR_ID == cnvr.id],
                          id = cnvr.id)
  cnvr.size <- cnvr.dat$end - cnvr.dat$start
  
  plot.label <- paste0( "chr", cnvr.dat$chr, ":", cnvr.dat$start, "-", cnvr.dat$end)
  
  # add a check to make sure  chrs match
  gene.dat <- annot.dat[ annot.dat$ID == cnvr.dat$id, ]
  gene.dat <- gene.dat[ !duplicated(gene.dat$name2), ]
  
  if( cnvr.dat$chr != unique(gene.dat$Chr) )
  {
    print(paste0("cnvr.dat is on chr", cnvr.dat$chr, " but gene data is mapped to chr", gene.dat$Chr[1], "."))
    stop("are you using the correct CNVR_ID?")
  }
  
  
  df <- data.frame( start = gene.dat$Start, end = gene.dat$End, name = gene.dat$name2)
  df$position <- seq(1:length(df$start)) + 1
  
  xmin <- min(c(cnvr.dat$start, df$start))  - (cnvr.size/4)
  xmax <- max(c(cnvr.dat$end, df$end)) + 10000
  
  p1 <- ggplot(data = df) + 
    geom_segment(aes(x = start, xend = end, y = position, yend = position), size = 4, color = "blue") +
    geom_segment(data = cnvr.dat, aes(x = start, xend = end, y = 1,  yend = 1), size = 4, color = "red") +
    geom_label(data = df, aes(x = start, y = position, label = name), size = 4, nudge_x = -nchar(df$name) * 0.02 * cnvr.size) +
    geom_label(data = cnvr.dat, aes(x  = start +  (end - start)/2, y = 1, label = id), nudge_y = max(df$position) * 0.025) +
    geom_vline(xintercept = cnvr.dat$start, linetype = "dashed", color = "red") +
    geom_vline(xintercept = cnvr.dat$end, linetype = "dashed", color = "red") +
    xlim( xmin, xmax) + 
    ggtitle(plot.label) +
    theme_minimal() 
  
  return(p1)
}
# ------------------------------------------------------------------------------ #
