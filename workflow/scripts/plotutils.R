# the utilties
library(ggplot2)
library(scales)
library(RColorBrewer)
# library(dplyr)
library(grid)
# library(gridBase)
library(gridExtra)
library(data.table)
library(gtable)
library(multidplyr)
library(dplyr)
library(tidyr)
# source("http://bioconductor.org/biocLite.R")
# biocLite("karyoploteR")
# BiocManager::install("karyoploteR")
library(karyoploteR)
library(GenomicRanges)
library(cowplot)
library(glue)
library(viridis)
library(stringr)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# BiocManager::install("HelloRanges")
library(HelloRanges)
# install.packages('scatterpie')
library(ggforce)
library(xlsx)
if (!require("ggnewscale")) install.packages("ggnewscale")


# chromosome names
CHRS <<- c(paste0("chr", seq(1, 22)), "chrX", "chrY", "chrM", "chrMT")
NOYM <- CHRS[!CHRS %in% c("chrY", "chrMT", "chrM")]
NOM <- CHRS[!CHRS %in% c("chrMT", "chrM")]

# colors to use
GRAY <- "#2F4F4F"
RED <- "#af0404"
BLUE <- "#3282b8"
NEWCOLOR <- RED
OLDCOLOR <- GRAY
COLORS <<- c(`T2T CHM13` = NEWCOLOR, GRCh37 = BLUE, GRCh38 = GRAY, `Celera WGSA` = "#ede682", `HG00733 pat` = "#96bb7c", `HG00733 mat` = "#ade498", `WGAC` = "#000000")
V <- "chm13.draft_v1.0_plus38Y"
ACHRO <<- paste0("chr", c(13, 14, 15, 21, 22))



readbed <- function(f, tag, rm = F, chrfilt = FALSE, fai = FAI) {
  df <- fread(glue(f))
  df
  colnames(df)[1:3] <- c("chr", "start", "end")
  df <- merge(df, fai[, c("chr", "chrlen")], by = "chr", all.x = TRUE)
  if (rm) {
    colnames(df) <- c("chr", "start", "end", "t1", "len", "strand", "type", "subtype", "x", "y")[1:length(colnames(df))]
  }
  if (chrfilt) {
    df <- df[chr %in% NOM & chr2 %in% NOM]
  }

  if ("chr2" %in% colnames(df)) {
    df$intra <- df$chr == df$chr2
    df$chr2 <- factor(df$chr2, levels = c(CHRS, unique(df$chr2[which(!df$chr2 %in% CHRS)])), ordered = TRUE)
    df <- df[chr %in% NOM & chr2 %in% NOM]
  }
  df$chr <- factor(df$chr, levels = c(CHRS, unique(df$chr[which(!df$chr %in% CHRS)])), ordered = TRUE)
  df$Assembly <- tag
  df <- df[order(chr, start)]
  df[chr != "chrM"]
}

grtodf <- function(gr) {
  df <- data.table(data.frame(gr))
  colnames(df)[1:3] <- c("chr", "start", "end")
  df$chr <- factor(df$chr, levels = c(CHRS, unique(df$chr[which(!df$chr %in% CHRS)])), ordered = TRUE)
  return(df)
}



#
# util funcations
#
add_genes <- function(df, all = FALSE, trim = FALSE, v1.1 = FALSE, genelist = GENES) {
  # df=rdg[,c(4,5,6,1)]
  # x = do_bedtools_intersect(toGRanges(df[,1:4]), toGRanges(GENES[,1:4]), loj = T); nrow(df)
  tmp <- genelist
  tmpall <- ALL_GENES
  if (v1.1) {
    tmp <- GENES_V1.1
    tmpall <- ALL_GENES_V1.1
  }

  if (all) {
    # fine ones with ORF
    ORF <- do.call(paste0, tmpall[, 1:12]) %in% do.call(paste0, tmp[, 1:12])
    sum(ORF)
    tmp <- tmpall
    tmp$ORF <- ORF
  } else {
    tmp$ORF <- T
  }

  o <- GenomicRanges::findOverlaps(toGRanges(df), toGRanges(tmp))
  if (trim) {
    tmp <- tmp[, c("gene", "ORF")]
  }
  cbind(df[queryHits(o)], tmp[subjectHits(o)])
}

merged_overlaps <- function(df, toadd) {
  toadd <- as.data.table(toadd)
  df <- as.data.table(df)
  o <- GenomicRanges::findOverlaps(
    toGRanges(df),
    GenomicRanges::reduce(toGRanges(toadd))
  )
  length(unique(queryHits(o)))
  length(queryHits(o))
  # nrow(df)
  cbind(df[queryHits(o)], toadd[subjectHits(o)])
}


has_genes <- function(df, all = FALSE) {
  # df=rdg[,c(4,5,6,1)]
  # x = do_bedtools_intersect(toGRanges(df[,1:4]), toGRanges(GENES[,1:4]), loj = T); nrow(df)
  tmp <- GENES
  if (all) {
    # fine ones with ORF
    ORF <- do.call(paste0, ALL_GENES[, 1:12]) %in% do.call(paste0, GENES[, 1:12])
    sum(ORF)
    tmp <- ALL_GENES
    tmp$ORF <- ORF
  } else {
    tmp$ORF <- Toverlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
    hits <- hits[percentOverlap > 0.5]
  }
  o <- GenomicRanges::findOverlaps(toGRanges(df), toGRanges(tmp))
  unique(queryHits(o))
}

is_sd <- function(df) {
  tmp <- GenomicRanges::reduce(toGRanges(SEDEF))
  o <- GenomicRanges::findOverlaps(toGRanges(df), tmp)
  unique(queryHits(o))
}


rgntag <- function(q, r, tag, minoverlap = 0, mincov = 0.5, rreduce = FALSE) {
  q <- as.data.table(q)
  r <- as.data.table(r)
  if (rreduce) {
    q <- data.table(data.frame(GenomicRanges::reduce(toGRanges(q))))
    colnames(q)[1:3] <- c("chr", "start", "end")
    q$chr <- factor(q$chr, levels = c(CHRS, unique(q$chr[which(!q$chr %in% CHRS)])), ordered = TRUE)
  }
  c <- do_bedtools_coverage(toGRanges(q[, 1:3]), toGRanges(r[, 1:3]))
  overlaps <- (c$fraction > mincov) & (c$covered > minoverlap)

  q[[tag]] <- overlaps
  q[[paste0("fraction_", tag)]] <- c$fraction
  q[[paste0("covered_", tag)]] <- c$covered
  q
}

overlaps <- function(q, r, minoverlap = 0, mincov = 0.0, rreduce = FALSE) {
  q <- as.data.table(q)
  r <- as.data.table(r)
  colnames(q)[1:3] <- c("chr", "start", "end")
  if (rreduce) {
    q <- data.table(data.frame(GenomicRanges::reduce(toGRanges(q))))
    q$chr <- factor(q$chr, levels = c(CHRS, unique(q$chr[which(!q$chr %in% CHRS)])), ordered = TRUE)
  }
  c <- do_bedtools_coverage(toGRanges(q[, 1:3]), toGRanges(r[, 1:3]))
  # print(c)
  overlaps <- (c$fraction > mincov) & (c$covered > minoverlap)
  return(overlaps)
}

overlap_either <- function(q, r, minoverlap = 250, mincov = 0.25, rreduce = FALSE) {
  # q=sedef; r=NEW
  o1 <- overlaps(q[, c("chr", "start", "end")], r, minoverlap = minoverlap, mincov = mincov, rreduce = FALSE)
  o2 <- overlaps(q[, c("chr2", "start2", "end2")], r, minoverlap = minoverlap, mincov = mincov, rreduce = FALSE)
  return(o1 | o2)
}

achro <- function(df) {
  cens <-
    df[(chr == "chr13" & end < CENS[chr == "chr13"]$start) |
      (chr == "chr14" & end < CENS[chr == "chr14"]$start) |
      (chr == "chr15" & start < CENS[chr == "chr15"]$start) |
      (chr == "chr21" & start < CENS[chr == "chr21"]$start) |
      (chr == "chr22" & start < CENS[chr == "chr22"]$start)]
}
is_achro <- function(df) {
  (df$chr == "chr13" & df$end < CENS[chr == "chr13"]$start) |
    (df$chr == "chr14" & df$end < CENS[chr == "chr14"]$start) |
    (df$chr == "chr15" & df$start < CENS[chr == "chr15"]$start) |
    (df$chr == "chr21" & df$start < CENS[chr == "chr21"]$start) |
    (df$chr == "chr22" & df$start < CENS[chr == "chr22"]$start)
}

chrplot <- function(df) {
  p <- ggplot(data = df[chr != "chrMT"], aes(x = start / 1000000, weight = (end - start))) +
    geom_density(adjust = 0.5) +
    scale_x_continuous(labels = comma) +
    facet_wrap(chr ~ ., ncol = 2) +
    xlab("Genomic position (Mbp)") +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    theme(strip.text = element_text(size = 10, margin = margin()), strip.background = element_blank())
  return(p)
}

addchr <- function(yy = 0) {
  return(geom_segment(data = fai, aes(x = 0, xend = start, y = yy, yend = yy)))
}

bedlength <- function(df, merge = TRUE) {
  gr <- GRanges(as.data.table(df))
  if (merge) {
    x <- sum(width(GenomicRanges::reduce(gr)))
  } else {
    x <- sum(width(gr))
  }
  x
}

find_num_sds <- function(df) {
  length(unique(df$unique_id))
}

gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


dupmasker_pal <- function() {
  d <- fread(glue("data/{V}_dupmasker_colors.bed"))
  d$hex <- sapply(strsplit(d$V9, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
  u <- unique(d[, c("V4", "hex")])
  y <- unlist(list(u$hex))
  names(y) <- u$V4
  y
}

in_rgn <- function(df, rgn) {
  # df = tbl_nof
  # rgn = "chr6:28,510,120-33,480,577"
  rgn <- gsub(",", "", rgn)
  x <- unlist(strsplit(rgn, "[:-]"))
  chr <- x[1]
  start <- as.integer(x[2])
  end <- as.integer(x[3])
  cond <- (chr == df$chr) & (start <= df$start & df$end <= end)
  # print(sum(cond))
  return(cond)
}


pairwise_plot <- function(df, facet = FALSE, inter = TRUE) {
  if (facet) {
    df <- df[(df$chr %in% unique(chm13$chr))]
  }
  if (!inter) {
    df <- df[df$chr %in% NOM & df$chr2 %in% NOM]
    df <- df[(df$chr == df$chr2)]
  }
  # [df=multiple_sds]
  df$cut <- cut(100 * df$fracMatch, breaks = seq(85.5, 100.5, 0.5), include.lowest = TRUE, right = FALSE)

  p1 <- ggplot(data = df) +
    geom_histogram(aes(100 * fracMatch, weight = alnB / 1e6, fill = Assembly),
      alpha = 0.9, color = "black", breaks = seq(89.5, 100.5, 0.5), position = "dodge", closed = "left"
    ) +
    theme_cowplot() +
    scale_y_continuous(labels = comma) +
    scale_x_continuous(breaks = seq(90, 100, 0.5), labels = c(90, "", 91, "", 92, "", 93, "", 94, "", 95, "", 96, "", 97, "", 98, "", 99, "", 100)) +
    ylab("Sum of aligned bases (Mbp)") +
    xlab("% identity of pairwise alignments") +
    theme(legend.position = "none") +
    guides(fill = guide_legend(ncol = length(COLORS))) +
    scale_fill_manual(values = COLORS)
  # theme(axis.text.x = element_text(size=7));p1



  p2 <- ggplot(data = df) +
    geom_histogram(aes(alnB, weight = alnB / 1e6, fill = Assembly),
      alpha = 0.85, color = "black", bins = 20, position = "dodge"
    ) +
    ylab("Sum of aligned bases (Mbp)") +
    xlab("Length of pairwise alignments") +
    scale_y_continuous(labels = comma) +
    scale_x_log10(labels = comma) +
    annotation_logticks(sides = "b") +
    theme_cowplot() +
    scale_fill_manual(values = COLORS) +
    theme(legend.position = "none")

  if (facet) {
    p1 <- p1 + facet_wrap(chr ~ ., ncol = 3)
    # p2 = p2 + facet_grid(chr ~ .) + theme(legend.position = "top")
    # p=plot_grid(p1, p2, ncol = 2,labels=c("a","b"))
    p <- plot_grid(p1, ncol = 1, labels = c("a"))
    return(p)
  }

  # p3 = ggplot(data = df) + geom_hex(aes(x=length, y=alnB/length)) +
  #  scale_x_log10() +facet_grid( assembly~.)

  # p=plot_grid(p1,p2,nrow = 2,labels=c("a","b"))
  print("here")
  p <- c()
  p$identity <- p1
  p$length <- p2
  return(p)
}

make_windows <- function(fai, window) {
  tmp <- data.table(fai %>% group_by(chr, chrlen) %>% summarise(start = list(seq(0, chrlen, window))) %>% unnest(col = start))
  tmp$end <- tmp$start + window
  tmp$end[tmp$end > tmp$chrlen] <- tmp$chrlen[tmp$end > tmp$chrlen] - 1
  tmp$chr <- factor(tmp$chr, levels = c(CHRS, unique(tmp$chr[which(!tmp$chr %in% CHRS)])), ordered = TRUE)
  return(tmp[, c("chr", "start", "end", "chrlen")])
}

zoom <- function(df, lift) {
  df$mod <- NA
  df$zname <- NA
  # for(x in lift$chr){
  for (i in 1:nrow(lift)) {
    row <- lift[i, ]
    cond <- overlaps(df, row)
    df$mod[cond] <- row$start
    df$zname[cond] <- row$name
  }
  df$start <- df$start - df$mod
  df$end <- df$end - df$mod
  if ("chr2" %in% colnames(df)) {
    df$mod2 <- NA
    df$zname2 <- NA
    for (x in lift$chr) {
      row <- lift[chr == x]
      cond <- overlaps(df[, c("chr2", "start2", "end2")], row)
      df$mod2[cond] <- row$start
      df$zname2[cond] <- row$name
    }
    df$start2 <- df$start2 - df$mod2
    df$end2 <- df$end2 - df$mod2
    df <- df[!is.na(start2)]
    df$start2[df$start2 < 0] <- 0
  }
  df <- df[!is.na(start)]
  df$start[df$start < 0] <- 0
  return(df)
}


length_stats <- function(df) {
  lengths <- rev(sort(df$length))
  total <- sum(lengths)
  cumm <- 0
  N50 <- 0
  for (len in lengths) {
    if (cumm >= total / 2) {
      N50 <- len
      break
    }
    cumm <- cumm + len
  }
  data.table(
    label = c("Count", "Total bp", "N50", "median", "mean"),
    value = c(length(lengths), total, N50, median(lengths), mean(lengths))
  )
}