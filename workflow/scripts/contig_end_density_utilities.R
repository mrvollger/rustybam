#!/usr/bin/env Rscript
if (!require("tidyverse"))
  install.packages("tidyverse")
if (!require("ggnewscale"))
  install.packages("ggnewscale")
if (!require("ggrepel"))
  install.packages("ggrepel")
if (!require("data.table"))
  install.packages("data.table")
if (!require("glue"))
  install.packages("glue")
if (!require("RColorBrewer"))
  install.packages("RColorBrewer")
if (!require("scales"))
  install.packages("scales")
if (!require("cowplot"))
  install.packages("cowplot")
if (!require("argparse"))
  install.packages("argparse")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!require("karyoploteR"))
  BiocManager::install("karyoploteR")
if (!require("GenomicRanges"))
  BiocManager::install("GenomicRanges")
if (!require("HelloRanges"))
  BiocManager::install("HelloRanges")
library(gridExtra)
library(openxlsx)
library(ggplotify)
library(ggridges)
library(ggpmisc)
library(ggforce)
library(ggpubr)

odir <<- "../../plots"
odir <<- "~/mvollger@uw.edu - Google Drive/Shared drives/AssemblyGaps/Figures/mrv_figures/."


my_ggsave <- function(filename, plot = last_plot(), ...) {
  filename = glue(filename)
  ggsave(filename, plot = plot, ...)
  basename = tools::file_path_sans_ext(filename)
  table_out = paste0(basename, ".datatable.tsv")
  data = apply(plot$data, 2, as.character)
  print(head(data))
  write.table(
    data,
    file = table_out,
    sep = "\t",
    row.names = F,
    quote = F
  )
}

#
# creat R data for the project to load and unload
#
if (F) {
  load("~/Desktop/Rdata/plotutils.data")
  fai_file <-
    "/Users/mrvollger/Desktop/EichlerVolumes/chm13_t2t/nobackups/assemblies/chm13_v1.1_plus38Y.fasta.fai"
  FAI <<-
    fread(fai_file, col.names = c("chr", "chrlen", "x", "y", "z"))
  FAI$chr <- factor(FAI$chr,
                    levels = c(CHRS, unique(FAI$chr[which(!FAI$chr %in% CHRS)])),
                    ordered = TRUE)
  
  pre <-
    "~/Desktop/EichlerVolumes/assembly_breaks/nobackups/rustybam_2021-08-16/reference_alignment/CHM13_V1.1"
  pre <-
    "~/Desktop/EichlerVolumes/assembly_breaks/nobackups/rustybam_2021-09-23/reference_alignment/CHM13_V1.1"
  inversions = readbed(
    glue(
      "{pre}/../../variants_freeze4inv_sv_inv_chm13_processed_arbigent_filtered_PAVgenAdded_onlyINVs_v1.0tov1.1.bed"
    ),
    tag = "INV"
  )
  inversion2 = fread(
    glue(
      "{pre}/../../variants_freeze4inv_sv_inv_chm13_processed_arbigent_filtered_PAVgenAdded_onlyINVs_v1.0tov1.1_mCNVoverlap.tsv"
    )
  )
  
  # header style for output tables
  hs <- createStyle(
    textDecoration = "BOLD",
    fontColour = "#000000",
    fontSize = 12,
    border = "bottom",
    borderStyle = "medium",
    halign = "center"
  )
  
  ### load seq content
  seq_content <-
    readbed(glue("{pre}/ends/all.ends.nuc.content.bed"),
            tag = "nuc")
  colnames(seq_content) <- gsub("\\d+_", "", colnames(seq_content))
  seq_content <- seq_content %>%
    mutate(pct_ga = (num_G + num_A) / seq_len,
           pct_tc = (num_T + num_C) / seq_len) %>%
    data.table()
  
  nuc <- readbed(glue("{pre}/ends/all.nuc.content.bed"),
                 tag = "nuc")
  ################################################
  
  link <-
    "https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/main/sample_metadata/hprc_year1_sample_metadata.txt"
  pop <- fread(link)
  
  read_data <-
    fread(
      "https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/main/read_metadata/hprc_year1_sample_level_misc_metadata_HiFi.csv"
    ) %>%
    rowwise() %>%
    mutate(
      G = 6.3e9,
      lambda = total_bp / G,
      expected_num_breaks = 2 * sum(dpois(seq(0, 4), lambda)) * G / mean
    ) %>%
    data.table()
  
  # read in the ends
  df <- fread(glue("{pre}/ends/all.ends.bed")) %>%
    drop_na() %>%
    separate(`reference_name`,
             into = c("sample", "hap", "tig"),
             sep = "#") %>%
    drop_na() %>%
    merge(pop,
          by.x = "sample",
          by.y = "Sample",
          all.x = T) %>%
    relocate("sample", .after = last_col()) %>%
    arrange(Superpopulation, sample, hap) %>%
    mutate(sample = factor(sample, levels = unique(sample)),
           ID = paste(sample, hap, tig)) %>%
    data.table()
  df$NID <- seq(nrow(df))
  
  df[is.na(Superpopulation)]$Superpopulation <- ""
  df[Superpopulation == ""]$Superpopulation <- "EUR"
  
  gdf <-
    toGRanges(df[, c("#query_name", "query_start", "query_end", "NID")])
  nonr_sd <- GenomicRanges::reduce(toGRanges(SEDEF_V1.1))
  alpha <-
    GenomicRanges::reduce(toGRanges(RM_V1.1[grep("Alpha", RM_V1.1$t1)]))
  othersat <-
    GenomicRanges::reduce(toGRanges(RM_V1.1[!grepl("Alpha", RM_V1.1$t1) &
                                              RM_V1.1$type == "Satellite"]))
  lowcom <-
    GenomicRanges::reduce(toGRanges(RM_V1.1[grep("Simple_repeat|Low_complexity", RM_V1.1$type)][, 1:3]))
  lowcom <- lowcom[width(lowcom) > 100]
  bedlength(lowcom) / 1e6
  
  #### add in the gc tontent
  o <- findOverlaps(gdf + 1e4, toGRanges(seq_content))
  bigdf <- cbind(df[queryHits(o)], seq_content[subjectHits(o)])
  df <- bigdf %>%
    group_by_at(vars(colnames(df))) %>%
    summarise(
      Max_GA_Frac = max(pct_ga),
      Max_TC_Frac = max(pct_tc),
      Max_GC_Frac = max(pct_gc),
      Max_AT_Frac = max(pct_at)
    ) %>%
    arrange(NID) %>%
    data.table()
  
  df$acro_p <-
    is_achro(data.table(
      chr = df$`#query_name`,
      start = df$query_start,
      end = df$query_end
    ))
  
  # add in the number of breaks in the regions to deinfe the liekely possion breaks
  # gdf <- toGRanges(df[, c("#query_name", "query_start", "query_end", "NID")])
  df$n_contig_ends_in_region <- countOverlaps(gdf + 1e5, gdf)
  sum(df$n_contig_ends_in_region <= 2)
  
  ### add in the breaks content
  df$sequence_context <- "other"
  df[(query_length - query_end < 2e5) |
       query_start < 2e5]$sequence_context <- "Chromosome end"
  df[n_contig_ends_in_region <= 2 &
       sequence_context == "other"]$sequence_context <-
    "Poisson breaks"
  
  df[overlaps(gdf, nonr_sd, mincov = 0.10) &
       (Max_GA_Frac >= 0.80 | Max_TC_Frac >= 0.80) &
       (sequence_context == "other")]$sequence_context <-
    "SD and High GA/TC (80%)"
  df[(Max_GA_Frac >= 0.80 |
        Max_TC_Frac >= 0.80) &
       (sequence_context == "other")]$sequence_context <-
    "High GA/TC (80%)"
  df[overlaps(gdf, nonr_sd, mincov = 0.10) &
       sequence_context == "other"]$sequence_context <- "SD"
  
  df[overlaps(gdf, alpha, mincov = 0.10) &
       sequence_context == "other"]$sequence_context <- "Alpha"
  df[overlaps(gdf, othersat, mincov = 0.10) &
       sequence_context == "other"]$sequence_context <- "Satellite"
  df[overlaps(gdf, lowcom, mincov = 0.10) &
       sequence_context == "other"]$sequence_context <-
    "10% Low Complexity"
  
  df[Max_GC_Frac >= 0.75 &
       (sequence_context == "other")]$sequence_context <-
    "High GC (75%)"
  df[Max_AT_Frac >= 0.80 &
       (sequence_context == "other")]$sequence_context <-
    "High AT (80%)"
  # df[`T2T HiFi low coverage` & (sequence_context == "other")]$sequence_context ="T2T HiFi low coverage"
  
  
  mycolors <- c(
    `Chromosome end` = "green",
    `Poisson breaks` = "cyan",
    `SD and High GA/TC (80%)` = "orange",
    `High GA/TC (80%)` = "darkgreen",
    SD = NEWCOLOR,
    Alpha = "purple",
    Satellite = "blue",
    `10% Low Complexity` = "lightblue",
    `High GC (75%)` = "green",
    `High AT (80%)` = "lightgreen",
    `T2T HiFi low coverage` = "orange",
    other = OLDCOLOR
  )
  df$sequence_context <-
    factor(df$sequence_context, rev(names(mycolors)))
  
  dfr <- fread(glue("{pre}/ends/windowed.all.ends.bed"))
  
  save.image(file = glue("{odir}/contig_ends.RData"))
} else {
  load(glue("{odir}/contig_ends.RData"))
}



make_karyoplot <-
  function(d1,
           d2,
           c1 = NEWCOLOR,
           c2 = OLDCOLOR,
           ym = 150,
           window.size = 1e6,
           chromosomes = NOM) {
    d1 <<- d1
    d2 <<- d2
    c1 <<- c1
    c2 <<- c2
    ym <<- ym
    chromosomes <<- chromosomes
    window.size <<- window.size
    pp <<- getDefaultPlotParams(plot.type = 2)
    pp$ideogramheight <<- 300
    ideo <- as.ggplot(
      expression(
        kp <- plotKaryotype(
          genome = GENOME,
          cytobands = CYTOFULL,
          chromosomes = chromosomes,
          plot.params = pp,
          plot.type = 2
        ),
        kpPlotRegions(
          kp,
          data = lowcom,
          data.panel = "ideogram",
          col = transparent("lightblue", .0),
          border = NA
        ),
        kpPlotRegions(
          kp,
          data = othersat,
          data.panel = "ideogram",
          col = transparent("blue", .0),
          border = NA
        ),
        kpPlotRegions(
          kp,
          data = alpha,
          data.panel = "ideogram",
          col = transparent("purple", .0),
          border = NA
        ),
        kpPlotRegions(
          kp,
          data = nonr_sd,
          data.panel = "ideogram",
          col = transparent(NEWCOLOR, .0),
          border = NA
        ),
        kpPlotDensity(
          kp,
          data = d1,
          data.panel = 1,
          ymax = ym,
          window.size = window.size,
          col = transparent(c1, 0.25)
        ),
        kpPlotDensity(
          kp,
          data = d2,
          col = transparent(c2, 0.25),
          window.size = window.size,
          ymax = ym * length(d1) / length(d2),
          data.panel = 2
        )
      )
    )
    ideo
  }

mycount <- function(vec) {
  # vec = sample_n(df,100)$sequence_context
  tbl <- sort(table(vec)[table(vec) >= 1], decreasing = T)
  paste(names(tbl), tbl, sep = ":", collapse = "; ")
}

add_contig_ends <- function(q.df, end.df, slop = 0) {
  go <- findOverlaps(toGRanges(q.df[, 1:3]) + slop, toGRanges(end.df))
  full <- cbind(q.df[queryHits(go)], end.df[subjectHits(go)])
  out <- full %>%
    filter(!(sequence_context %in% c("Chromosome end", "Poisson breaks"))) %>%
    group_by_at(vars(colnames(q.df))) %>%
    summarise(
      `# contig ends` = n(),
      n_samples = length(unique(sample)),
      n_haplotypes = length(unique(paste(sample, hap))),
      sequence_context = mycount(sequence_context)
    ) %>%
    ungroup() %>%
    arrange(-`# contig ends`) %>%
    filter(n_samples > 2) %>%
    data.table()
  out$acro_p <- is_achro(out)
  out
}
