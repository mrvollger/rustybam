#!/usr/bin/env Rscript
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggnewscale")) install.packages("ggnewscale")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("data.table")) install.packages("data.table")
if (!require("glue")) install.packages("glue")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("scales")) install.packages("scales")
if (!require("cowplot")) install.packages("cowplot")
if (!require("argparse")) install.packages("argparse")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!require("karyoploteR")) BiocManager::install("karyoploteR")
if (!require("GenomicRanges")) BiocManager::install("GenomicRanges")
if (!require("HelloRanges")) BiocManager::install("HelloRanges")

library(openxlsx)
library(ggplotify)
library(ggridges)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
odir <<- "../../plots"
my_ggsave <- function(filename, plot = last_plot(), ...){
  filename=glue(filename)
  ggsave(filename, plot=plot, ...)
  basename = tools::file_path_sans_ext(filename)
  table_out = paste0(basename, ".datatable.tsv")
  data=apply(plot$data,2,as.character)
  print(head(data))
  write.table(data, file=table_out, sep="\t", row.names = F, quote = F)
}

if (F) {
  load("~/Desktop/Rdata/plotutils.data")
  fai_file <- "/Users/mrvollger/Desktop/EichlerVolumes/chm13_t2t/nobackups/assemblies/chm13_v1.1_plus38Y.fasta.fai"
  FAI <<- fread(fai_file, col.names = c("chr", "chrlen", "x", "y", "z"))
  FAI$chr <- factor(
    FAI$chr,
    levels = c(CHRS, unique(FAI$chr[which(!FAI$chr %in% CHRS)])),
    ordered = TRUE
  )

  pre <- "~/Desktop/EichlerVolumes/assembly_breaks/nobackups/rustybam_2021-08-16/reference_alignment/CHM13_V1.1"
  pre <- "~/Desktop/EichlerVolumes/assembly_breaks/nobackups/rustybam_2021-09-23/reference_alignment/CHM13_V1.1"
  
  hs <- createStyle(
    textDecoration = "BOLD", fontColour = "#000000", fontSize = 12,
    border = "bottom", borderStyle = "medium", halign = "center"
  )

  # hifi_issues <- readbed(
  #  "https://raw.githubusercontent.com/marbl/CHM13-issues/main/v1.1/hifi.pri.issues.bed",
  #  "HiFi sequencing issues"
  # )
  # colnames(hifi_issues)[4] <- "HiFi notes"
  # hifi_issues

  ### load seq content
  seq_content <- readbed(
    glue("{pre}/ends/all.ends.nuc.content.bed.gz"),
    tag = "nuc"
  )
  colnames(seq_content) <- gsub("\\d+_", "", colnames(seq_content))
  seq_content <- seq_content %>%
    mutate(
      pct_ga = (num_G + num_A) / seq_len,
      pct_tc = (num_T + num_C) / seq_len
    ) %>%
    data.table()
 
  nuc <- readbed(
    glue("{pre}/ends/all.nuc.content.bed.gz"),
    tag = "nuc"
  )
    ################################################

  
  link <- "https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/main/sample_metadata/hprc_year1_sample_metadata.txt"
  pop <- fread(link)

  read_data <- fread("https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/main/read_metadata/hprc_year1_sample_level_misc_metadata_HiFi.csv") %>%
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
    separate(`reference_name`, into = c("sample", "hap", "tig"), sep = "#") %>%
    drop_na() %>%
    merge(pop, by.x = "sample", by.y = "Sample", all.x = T) %>%
    relocate("sample", .after = last_col()) %>%
    arrange(Superpopulation, sample, hap) %>%
    mutate(
      sample = factor(sample, levels = unique(sample)),
      ID = paste(sample, hap, tig)
    ) %>%
    data.table()
  df$NID <- seq(nrow(df))

  df[is.na(Superpopulation)]$Superpopulation <- ""
  df[Superpopulation == ""]$Superpopulation <- "EUR"

  gdf <- toGRanges(df[, c("#query_name", "query_start", "query_end", "NID")])
  nonr_sd <- GenomicRanges::reduce(toGRanges(SEDEF_V1.1))
  alpha <- GenomicRanges::reduce(toGRanges(RM_V1.1[grep("Alpha", RM_V1.1$t1)]))
  othersat <- GenomicRanges::reduce(toGRanges(RM_V1.1[!grepl("Alpha", RM_V1.1$t1) & RM_V1.1$type == "Satellite"]))
  lowcom <- GenomicRanges::reduce(toGRanges(RM_V1.1[grep("Simple_repeat|Low_complexity", RM_V1.1$type)][, 1:3]))
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

  df$acro_p <- is_achro(data.table(chr = df$`#query_name`, start = df$query_start, end = df$query_end))

  # add in the number of breaks in the regions to deinfe the liekely possion breaks
  df$n_contig_ends_in_region <- countOverlaps(gdf + 1e5, gdf)
  sum(df$n_contig_ends_in_region <= 2)

  ### add in the breaks content
  df$sequence_context <- "other"
  df[(query_length - query_end < 2e5) | query_start < 2e5]$sequence_context <- "Chromosome end"
  df[n_contig_ends_in_region <= 2 & sequence_context == "other"]$sequence_context <- "Poisson breaks"

  df[overlaps(gdf, nonr_sd, mincov = 0.10) &
    (Max_GA_Frac >= 0.80 | Max_TC_Frac >= 0.80) &
    (sequence_context == "other")]$sequence_context <- "SD and High GA/TC (80%)"
  df[(Max_GA_Frac >= 0.80 | Max_TC_Frac >= 0.80) & (sequence_context == "other")]$sequence_context <- "High GA/TC (80%)"
  df[overlaps(gdf, nonr_sd, mincov = 0.10) & sequence_context == "other"]$sequence_context <- "SD"

  df[overlaps(gdf, alpha, mincov = 0.10) & sequence_context == "other"]$sequence_context <- "Alpha"
  df[overlaps(gdf, othersat, mincov = 0.10) & sequence_context == "other"]$sequence_context <- "Satellite"
  df[overlaps(gdf, lowcom, mincov = 0.10) & sequence_context == "other"]$sequence_context <- "10% Low Complexity"

  df[Max_GC_Frac >= 0.75 & (sequence_context == "other")]$sequence_context <- "High GC (75%)"
  df[Max_AT_Frac >= 0.80 & (sequence_context == "other")]$sequence_context <- "High AT (80%)"
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
  df$sequence_context <- factor(df$sequence_context, rev(names(mycolors)))

  dfr <- fread(glue("{pre}/ends/windowed.all.ends.bed"))

  save.image(file = glue("{odir}/contig_ends.RData"))
} else {
  load(glue("{odir}/contig_ends.RData"))
}
########################################################################

make_karyoplot <- function(d1, d2, c1 = NEWCOLOR, c2 = OLDCOLOR, ym = 150, window.size = 1e6, chromosomes = NOM) {
  d1 <<- d1
  d2 <<- d2
  c1 <<- c1
  c2 <<- c2
  ym <<- ym
  chromosomes <<- chromosomes
  window.size <<- window.size
  pp <<- getDefaultPlotParams(plot.type = 2)
  pp$ideogramheight <<- 300
  ideo <- as.ggplot(expression(
    kp <- plotKaryotype(
      genome = GENOME, cytobands = CYTOFULL,
      chromosomes = chromosomes,
      plot.params = pp, plot.type = 2
    ),
    kpPlotRegions(kp,
      data = lowcom, data.panel = "ideogram",
      col = transparent("lightblue", .0), border = NA
    ),
    kpPlotRegions(kp,
      data = othersat, data.panel = "ideogram",
      col = transparent("blue", .0), border = NA
    ),
    kpPlotRegions(kp,
      data = alpha, data.panel = "ideogram",
      col = transparent("purple", .0), border = NA
    ),
    kpPlotRegions(kp,
      data = nonr_sd, data.panel = "ideogram",
      col = transparent(NEWCOLOR, .0), border = NA
    ),
    kpPlotDensity(kp,
      data = d1,
      data.panel = 1,
      ymax = ym,
      window.size = window.size,
      col = transparent(c1, 0.25)
    ),
    kpPlotDensity(kp,
      data = d2,
      col = transparent(c2, 0.25),
      window.size = window.size,
      ymax = ym * length(d1) / length(d2),
      data.panel = 2
    )
  ))
  ideo
}

exclude <- c("Chromosome end", "Poisson breaks", "High GA/TC (80%)", "SD and High GA/TC (80%)", "SD")
ideo <- make_karyoplot(
  gdf[df$sequence_context %in% c("SD and High GA/TC (80%)", "SD")],
  gdf[!(df$sequence_context %in% exclude)]
) +
  theme_map(font_size = 16) + ggtitle("SD contig ends vs non-SD")
ideo

df$simple_seq_content <- recode(df$sequence_context,
  `Chromosome end` = "Chromosome end or Poisson",
  `Poisson breaks` = "Chromosome end or Poisson",
  `SD and High GA/TC (80%)` = "SD and High GA/TC (80%)",
  `SD` = "SD and High GA/TC (80%)",
  `Alpha` = "Satellite",
  `Satellite` = "Satellite",
  `10% Low Complexity` = "Satellite",
  `High GC (75%)` = "other",
  `High AT (80%)` = "other",
  `other` = "other"
)
simple_colors <- c(
  `SD and High GA/TC (80%)` = NEWCOLOR,
  Satellite = "darkblue",
  `Chromosome end or Poisson` = "lightgreen",
  `High GA/TC (80%)` = "darkgreen",
  other = OLDCOLOR
)
df$simple_seq_content <- factor(df$simple_seq_content, levels = rev(names(simple_colors)))


simple_p <- ggplot(data = df, aes(y = simple_seq_content, fill = simple_seq_content)) +
  geom_bar(alpha = 0.9) +
  geom_text(stat = "count", aes(label = comma(..count..)), hjust = -0.5) +
  scale_fill_manual("", values = simple_colors) +
  ylab("") +
  xlab("Number of contig ends") +
  theme_cowplot() +
  theme(legend.position = "none") +
  theme(
    plot.margin = margin(t = 0.2, b = 0.2, l = 0.2, r = 2, "cm")
  ) +
  coord_cartesian(clip = "off") +
  ggtitle("Contig ends (gaps) in HPRC Hifiasm assemblies")
simple_p
my_ggsave("{odir}/simple_hifiasm_contig_ends_hist.pdf", height = 4, width = 8, plot = simple_p)


p <- ggplot(data = df, aes(y = sequence_context, fill = sequence_context)) +
  geom_bar(alpha = 0.9) +
  geom_text(stat = "count", aes(label = comma(..count..)), hjust = -0.5) +
  scale_fill_manual("", values = mycolors) +
  ylab("") +
  xlab("Number of contig ends") +
  theme_cowplot() +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off") + 
  theme(plot.margin=margin(r=1, unit="in")) +
  ggtitle("Contig ends (gaps) in HPRC Hifiasm assemblies")
p
my_ggsave(plot=p, height=9/2, width=16/2, file="{odir}/hifiasm_contig_ends_hist.pdf")

pp <- ggplot(data = df) +
  geom_bar(aes(y = paste(sample, hap), fill = Superpopulation)) +
  theme_cowplot() +
  theme(legend.position = "top") +
  scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
  ylab("") +
  xlab("Number of contig ends per sample")
pp
my_ggsave(plot=pp, height=12, width=8, file="{odir}/num_hprc_contig_ends.pdf")

z <- cowplot::plot_grid(p, pp, rel_heights = c(1, 3), ncol = 1)
z
f <- cowplot::plot_grid(z, ideo, ncol = 2, rel_widths = c(1, 2))
f
ggsave(glue("{odir}/1_hifiasm_contig_ends.pdf"), height = 16, width = 24, plot = f)


#-----------------------AFR-----------------------------------#

afr <- make_karyoplot(gdf[df$Superpopulation == "AFR"], gdf[df$Superpopulation != "AFR"], c1 = "orange", c2 = "blue") +
  theme_map(font_size = 16) + ggtitle("African contig end vs non-African")
my_ggsave("{odir}/4_afr_contig_ends.pdf", height = 12, width = 12)

#-----------------------TABLE-----------------------------------#


dfr <- dfr %>%
  filter((query_start > 1e5) & (query_length - query_end > 1e5)) %>%
  separate(`reference_name`, into = c("sample", "hap", "tig"), sep = "#") %>%
  group_by(`#chr`, start, end) %>%
  summarise(num_contig_breaks = n()) %>%
  # haps=list(hap)) %>%#, samples=list(samples), tigs=list(tig)) %>%
  ungroup() %>%
  mutate(percentile = ecdf(num_contig_breaks)(num_contig_breaks)) %>%
  data.table() %>%
  add_genes(trim = TRUE, genelist = GENES_V1.1[overlaps(GENES_V1.1, nonr_sd)]) %>%
  group_by(`#chr`, start, end, num_contig_breaks, percentile) %>%
  summarise(gene = paste0(unique(gene), collapse = ";")) %>%
  rgntag(nonr_sd, "SD", mincov = 0) %>%
  rgntag(c(alpha, othersat, lowcom), "Sat", mincov = 0) %>%
  relocate("gene", .after = last_col()) %>%
  data.frame() %>%
  arrange(-percentile) %>%
  data.table()
dfr <- dfr[, -c("SD", "Sat", "covered_Sat", "covered_SD")]
dfr

#-----------------------OTHER-----------------------------------#

o <- findOverlaps(toGRanges(RM_V1.1), gdf[df$sequence_context == "other"])
otherrm <- RM_V1.1[unique(queryHits(o))]
zz <- cowplot::plot_grid(ggplot(data = otherrm) +
  geom_bar(aes(y = type, weight = end - start)) +
  ggtitle("Annotations in other contig ends"),
ggplot(data = RM_V1.1) +
  geom_bar(aes(y = type, weight = end - start)) +
  ggtitle("Annotations whole genome"),
ncol = 1
)
other_ideo <- make_karyoplot(gdf[df$sequence_context == "other"], GRanges(), c1 = OLDCOLOR, window.size = 5e5, ym = 10) +
  ggtitle("Position of contig ends that were 'other'")
ggsave(file = glue("{odir}/3_annotations_in_other_contig_ends.pdf"), height = 12, width = 16, plot = cowplot::plot_grid(zz, other_ideo, ncol = 2))


ofemale <- gdf[!is.na(df$Sex) & df$sequence_context == "other" & df$Sex == "female"]
omale <- gdf[!is.na(df$Sex) & df$sequence_context == "other" & df$Sex == "male"]

make_karyoplot(ofemale, omale, ym = 40, window.size = 3e5, c1 = "orange", c2 = "blue", chromosomes = c("chrX")) +
  ggtitle("Female vs male contig ends that were 'other'")

nhaps <- table(unique(df[, c("sample", "Sex")])$Sex)

female <- !is.na(df$Sex) & df$sequence_context == "other" & df$Sex == "female" & df$`#query_name` == "chrX"
sum(female) / (nhaps["female"] * 2)
male <- !is.na(df$Sex) & df$sequence_context == "other" & df$Sex == "male" & df$`#query_name` == "chrX"
sum(male) / nhaps["male"]


#----------------------- genes gaps -----------------------------------#
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
      n_haplotypes = length(unique(paste(sample,hap))),
      sequence_context = mycount(sequence_context)
    ) %>%
    ungroup() %>%
    arrange(-`# contig ends`) %>%
    filter(n_samples > 2) %>%
    data.table()
  out$acro_p <- is_achro(out)
  out
}
simple_genes <- ALL_GENES_V1.1 %>%
  group_by(chr, gene) %>%
  summarise(start = min(start), end = max(end)) %>%
  filter(!grepl("\\d+\\.\\d+", gene)) %>%
  relocate(gene, .after = last_col()) %>%
  data.table()

gene_break_summary <- add_contig_ends(simple_genes, df, slop = 1e3)
gene_break_summary[grep("NPIPA9", gene_break_summary$gene)]
# openxlsx::write.xlsx(
#  list(`Contig ends within genes` =  gene_break_summary, `Genic windows with contig ends`=dfr),
#  headerStyle = hs,numFmt = "COMMA", colWidths="auto", gridLines=FALSE, colNames = TRUE,
#  "contig_ends.xlsx")

#-----------------------SD clusters -----------------------------------#
large_sd_blocks <- as.data.table(GenomicRanges::reduce(nonr_sd[width(nonr_sd) > 5e4] + 5e4))
colnames(large_sd_blocks)[1] <- "chr"
large_sd_blocks[start < 0]$start <- 0
large_sd_blocks <- add_genes(large_sd_blocks, genelist = simple_genes)[, c(1:3, 9)] %>%
  # drop the extra columns
  group_by(chr, start, end) %>%
  summarise(genes = list(sort(unique(gene)))) %>%
  data.table()
unique(df$sample)
large_sd_blocks_with_ends <- add_contig_ends(large_sd_blocks, df) %>%
  mutate(
    Mbp = (end - start) / 1e6,
    `contig ends per kbp` = `# contig ends` / (end - start) * 1000
  ) %>%
  arrange(chr, start) %>%
  arrange(-`# contig ends`) %>%
  relocate(genes, .after = last_col()) %>%
  data.table()


large_sd_blocks_with_ends
#save(file="~/Desktop/large_sd_blocks.Rdata", large_sd_blocks_with_ends)
#load(file="~/Desktop/large_sd_blocks.Rdata")
large_sd_blocks_with_ends$length = large_sd_blocks_with_ends$end - large_sd_blocks_with_ends$start
large_sd_blocks_with_ends$y=(large_sd_blocks_with_ends$length)
large_sd_blocks_with_ends$x=(large_sd_blocks_with_ends$`# contig ends`)
large_sd_blocks_with_ends$top_context = gsub(":.*", "", large_sd_blocks_with_ends$sequence_context)
large_sd_blocks_with_ends$genes=sapply(large_sd_blocks_with_ends$genes, paste, collapse=",")
dim(large_sd_blocks_with_ends)
library(ggpubr)
p.len_vs_n_samples = ggplot(data=large_sd_blocks_with_ends,
                            aes(x=n_haplotypes, y=y)
                            ) +
  geom_point()+
  scale_y_continuous(trans="log10", labels = comma)+
  #scale_x_continuous(trans="log10", labels = comma)+
  annotation_logticks(sides = "l")+
  geom_smooth(se=FALSE, method = "lm")+
  #facet_wrap(~top_context)+
  stat_cor( method="pearson")+
  theme_cowplot()+
  xlab("# broken assemblies in SD loci across the 94 HPRC haplotypes")+
  ylab("Length of the SD locus")
p.len_vs_n_samples  
my_ggsave(file="{odir}/len_vs_nsamples.pdf", plot=p.len_vs_n_samples, height=8, width=8)

#----------------------- windowed anlysis -----------------------------------#
myseqlengths <- FAI_v1.1$chrlen
names(myseqlengths) <- FAI_v1.1$chr
windows <- as.data.table(tileGenome(myseqlengths, tilewidth = 1e6, cut.last.tile.in.chrom = T))
colnames(windows)[1] <- "chr"
windows[start < 0]$start <- 0
windows <- add_genes(windows, genelist = simple_genes)[, c(1:3, 9)] %>%
  # drop the extra columns
  group_by(chr, start, end) %>%
  summarise(genes = list(sort(unique(gene)))) %>%
  data.table()
windows_with_ends <- add_contig_ends(windows, df) %>%
  mutate(
    Mbp = (end - start) / 1e6,
    `contig ends per kbp` = `# contig ends` / (end - start) * 1000
  ) %>%
  arrange(chr, start) %>%
  arrange(-`# contig ends`) %>%
  relocate(genes, .after = last_col()) %>%
  data.table()

#----------------------- Table -----------------------------------#
openxlsx::write.xlsx(
  list(
    `Contig ends in large SD blocks` = large_sd_blocks_with_ends,
    `Contig ends within genes` =  gene_break_summary,
    `1 Mbp windows with contig ends` = windows_with_ends
  ),
  headerStyle = hs, numFmt = "COMMA", colWidths = "auto", gridLines = FALSE, colNames = TRUE, overwrite = TRUE,
  glue("{odir}/contig_ends.xlsx")
)

#----------------------- Seq content genome wide -----------------------------------#
if (F) {
  if (T) {
    colnames(nuc) <- gsub("\\d+_", "", colnames(nuc))
    nuc$pct_ga <- (nuc$num_A + nuc$num_G) / nuc$seq_len
    nuc$pct_tc <- (nuc$num_T + nuc$num_C) / nuc$seq_len
  }
  AT <- NULL
  GA <- NULL
  GC <- NULL
  TC <- NULL
  for (i in seq(1e5)) {
    consecutive_windows_to_observe <- 20 * 10
    random_start <- runif(1, min = 1, max = nrow(nuc) - consecutive_windows_to_observe)
    rands <- nuc[random_start:(random_start + consecutive_windows_to_observe)]
    AT <- c(AT, max(rands$pct_at))
    GC <- c(GC, max(rands$pct_gc))
    GA <- c(GA, max(rands$pct_ga))
    TC <- c(TC, max(rands$pct_tc))
  }
  rand_seq <- data.table(AT, GC, GA, TC)
  rand_seq <- pivot_longer(rand_seq, cols = c("AT", "GC", "GA", "TC"))
  rand_seq %>%
    group_by(name) %>%
    summarise(sum(value > 0.8) / 1e5 * 100)

  rand_seq$name <- factor(rand_seq$name, levels = unique(rand_seq$name))
  p_seq_content <- ggplot(data = rand_seq) +
    geom_density_ridges(aes(x = value, fill = name, y = name), alpha = 0.5) +
    ylab("Sequence content") +
    xlab("Fraction of the 1,000 bp window") +
    scale_x_continuous(breaks = seq(0, 1, 0.05)) +
    theme_minimal_vgrid() +
    theme(legend.position = "none") +
    ggtitle("Simulation of sequnece content in 1,000 bp widnows (n=100,000)")
  p_seq_content
  my_ggsave("{odir}/2_simulation_of_sequence_content.pdf", plot = p_seq_content, height = 9, width = 12)
}


#----------------------- SD simulations wide -----------------------------------#
if (F) {
  high_nuc <- nuc[(num_G + num_A > 800) | (num_C + num_T > 800)]
  myseqlengths = FAI_v1.1$chrlen
  names(myseqlengths) = FAI_v1.1$chr
  windows <- tileGenome(myseqlengths, tilewidth = 1e3, cut.last.tile.in.chrom = T)
  smalldf <- df[!(sequence_context %in% c("Chromosome end", "Poisson breaks"))]
  n <- nrow(smalldf)
  sim <- NULL
  for (i in seq(10000)) {
    rand_idx <- runif(n, min = 1, max = length(windows))
    rand_windows <- windows[rand_idx]
    z <- findOverlaps(rand_windows, nonr_sd)
    s <- length(unique(queryHits(z)))

    zz <- findOverlaps(rand_windows, alpha)
    a <- length(unique(queryHits(zz)))

    zzz <- findOverlaps(rand_windows, othersat)
    o <- length(unique(queryHits(zzz)))

    zzzz <- findOverlaps(rand_windows, toGRanges(high_nuc) + 10e3)
    hn <- length(unique(queryHits(zzzz)))

    sim[[i]] <- data.table(
      name = c("SD", "Alpha", "Satellite", "High GA/TC (80%)"),
      value = c(s, a, o, hn)
    )
  }
  sim <- data.table::rbindlist(sim)
  sim_median <- sim %>%
    group_by(name) %>%
    summarise(median_simulation = median(value))

  obs <- data.table(
    name = c("SD", "Alpha", "Satellite", "High GA/TC (80%)"),
    value = c(
      sum(overlaps(smalldf, nonr_sd, mincov = 0.1)),
      sum(overlaps(smalldf, alpha, mincov = 0.1)),
      sum(overlaps(smalldf, othersat, mincov = 0.1)),
      sum(overlaps(smalldf, toGRanges(high_nuc) + 10e3))
    )
  ) %>% merge(sim_median) %>%
  mutate( observation = value, `fold increase`=value/median_simulation) %>%
  data.table()
  obs
  library(gridExtra)
  p_sd_content <- ggplot(data = sim) +
    geom_density_ridges(aes(x = value, y = name, fill = name), alpha = 0.5) +
    geom_point(data = obs, aes(x = value, y = name, color = name), size = 3) +
    geom_text(
      data = obs, aes(x = value, y = name, label = comma(value), color = name),
      vjust = -1
      ) +
    geom_text(data = obs, 
              aes(x = (value+median_simulation)/2, y = name, label= paste0("fold increase = ", round(`fold increase`, 2)), color = name),
              vjust=-3
              )+
    ylab("Sequence content") +
    xlab("# of windows with the annotation") +
    theme_minimal_grid() +
    theme(legend.position = "none") +
    scale_fill_manual(values = mycolors) +
    scale_color_manual(values = mycolors) +
    scale_x_continuous(labels = comma) +
    coord_cartesian(xlim = c(1, NA)) +
    annotation_custom(tableGrob(obs %>% dplyr::select(-value), 
                                theme=ttheme_minimal(), rows = NULL),
                      xmin = 10e3, xmax = 20e3, ymin = "Satellite", ymax = "Satellite") +
    ggtitle(glue("Simulation showing the number of windows (1 kbp width) out of {nrow(df)}\nthat overlap with seq class (n=10,000)"))
  p_sd_content
  obs
  my_ggsave("{odir}/2_simulation_of_SDs.pdf", plot = p_sd_content, height = 9, width = 12)
}


#----------------------- SD features of SD breaks -----------------------------------#
sd_breaks <- df[sequence_context == "SD", ]
sdo <- findOverlaps(toGRanges(sd_breaks), toGRanges(SEDEF_V1.1))
windows <- tileGenome(myseqlengths, tilewidth = 1e3, cut.last.tile.in.chrom = T)

extra_sd_breaks <- cbind(sd_breaks[queryHits(sdo)], SEDEF_V1.1[subjectHits(sdo)]) %>%
  group_by_at(vars(colnames(sd_breaks))) %>%
  summarise(
    fracMatchZ = max(fracMatch),
    max_frac_len = matchB[which.max(fracMatch)],
    matchB = max(matchB),
    fracMatch = fracMatch[which.max(matchB)]
  ) %>%
  data.table()
# extra_sd_breaks =  SEDEF_V1.1[unique(subjectHits(sdo))]#cbind(sd_breaks[queryHits(sdo)], SEDEF_V1.1[subjectHits(sdo)])
nrow(extra_sd_breaks)
sd.plot.df <- rbindlist(
  l = list(
    extra_sd_breaks[, c("fracMatch", "matchB")],
    SEDEF_V1.1[SEDEF_V1.1$original, c("fracMatch", "matchB")]
  ),
  idcol = "ContigEnd"
)
sd.plot.df$ContigEnd <- sd.plot.df$ContigEnd == 1

summary <- sd.plot.df %>%
  group_by(ContigEnd) %>%
  summarise(
    meanID = round(mean(fracMatch) * 100, 2),
    meanLen = comma(round(mean(matchB), 2)),
    medianID = round(median(fracMatch) * 100, 2),
    medianLen = comma(median(matchB))
  )

sd_col <- c(NEWCOLOR, OLDCOLOR)
names(sd_col) <- c(TRUE, FALSE)
sd_density <- ggplot(data = sd.plot.df, aes(x = matchB, y = fracMatch * 100, color = ContigEnd, fill = ContigEnd)) +
  annotation_custom(tableGrob(t(summary), theme=ttheme_minimal()), xmin = 6, xmax = 7, ymin = 91, ymax = 93) +
  geom_point(alpha = 0.25, size = 0.1) +
  geom_density_2d(size = 0.5, alpha = 0.75) +
  scale_fill_manual(values = sd_col) +
  scale_color_manual(values = sd_col) +
  xlab("Alignment length") +
  ylab("% identity") +
  scale_x_log10(labels = comma) +
  annotation_logticks(sides = "b") +
  ggtitle("Distribution of SD length and identity\nin contig ends vs whole genome") +
  theme_cowplot(font_size = 18)
my_ggsave("{odir}/5_SD_length_frac_density.pdf", height = 4 * 1.5, width = 6 * 1.5, plot = sd_density)
sd_density



ggplot() +
  geom_point(aes(x = c(0, 1), y = c(90, 100))) +
  annotation_custom(tableGrob(summary), xmin = 0, xmax = 1, ymin = 90, ymax = 93)






#----------------------- SD features of SD breaks -----------------------------------#

SMN <- large_sd_blocks_with_ends[grepl("SMN", large_sd_blocks_with_ends$genes), ]
LPA <- large_sd_blocks_with_ends[grepl("LPA", large_sd_blocks_with_ends$genes), ]

go <- findOverlaps(toGRanges(df), toGRanges(simple_genes))
with_genes <- cbind(df[queryHits(go)], simple_genes[subjectHits(go)])
add_in <- df$NID[!(df$NID %in% queryHits(go))]
with_genes <- rbindlist(list(with_genes, df[add_in]), fill = T)


plot.df <- with_genes[overlaps(with_genes, SMN)]

ggplot(data = plot.df) +
  geom_segment(aes(y = " genes", yend = " genes", x = start, xend = end, color = (gene == "SMN1" | gene == "SMN2"))) +
  geom_point(aes(x = query_start, y = Superpopulation)) +
  theme_cowplot() +
  scale_x_continuous(labels = comma) +
  theme(legend.position = "none")

df[, c(1:3, 7, 6, 5, 1:3)]



#----------------------- contig ends vs coverage -----------------------------------#
problems <- c(
  "SD and High GA/TC (80%)",
  "High GA/TC (80%)"
 )

ends_and_depht <- df %>%
  filter(sequence_context %in% problems) %>%
  group_by(sample, sequence_context) %>%
  summarise(Number_of_GA_TC_breaks = sum(sequence_context %in% problems), Sex = unique(Sex), Superpopulation = unique(Superpopulation)) %>%
  merge(read_data, by.x = "sample", by.y = "sample_id") %>%
  mutate(sequence_context = as.character(sequence_context)) %>%
  data.table() 
library(ggforce)
p.cov.breaks <- ggplot(
  data = ends_and_depht,
  aes(
    x = total_Gbp, 
    y = Number_of_GA_TC_breaks
    )
  ) +
  geom_point() +
  geom_smooth(method = "lm", alpha = 0.20, size = 0.25, se=FALSE) +
  geom_text_repel(aes(label=sample))+
  theme_cowplot() +
  scale_x_log10()+scale_y_log10()+annotation_logticks()+
  #facet_zoom(x = total_Gbp < 200, show.area=TRUE)+
  facet_row(~sequence_context, scales = "free")+
  ggtitle("Number of contig ends in high GA/TC (80%) vs total Gbp of sequencing") +
  theme(legend.position="bottom") + 
  xlab("Total Gbp of Hifi data") + ylab("# of GA/TC (80%) associated contig ends")
p.cov.breaks
my_ggsave("{odir}/7_sample_coverage_and_breaks.pdf", height = 6, width = 12, plot = p.cov.breaks)

read_data

p.read.n50 = ggplot(
  data = read_data[sample_id != "HG002"],
  aes(x = total_Gbp, y = N50)
  ) +
  geom_point() +
  geom_smooth(method = "lm", alpha = 0.20, size = 0.25, se=F) +
  theme_cowplot() +
  #facet_zoom(x = total_Gbp < 200, show.area=TRUE)+
  ggtitle("Read N50 vs total Gbp of sequencing")
p.read.n50
# facet_row(~Sex);p.cov.breaks
my_ggsave("{odir}/8_read_n50_vs_total_gbp.pdf", height = 9, width = 16, plot = p.read.n50)


#----------------------- collapsed bases -----------------------------------#
collapse <- fread("~/Desktop/EichlerVolumes/assembly_breaks/nobackups/hifi_read_alignments/results/collapse.bed")

col_sum <- collapse %>%
  group_by(sample) %>%
  summarise(`Collapsed bases` = sum(end - start), `Predicted number of unassembled bases` = sum((end - start) * coverage / median)) %>%
  pivot_longer(cols = c("Collapsed bases", "Predicted number of unassembled bases"))
sum2 <- col_sum %>%
  group_by(name) %>%
  summarise(mean = mean(value))

col.plot <- ggplot(data = col_sum, aes(y = sample, x = value, fill = name)) +
  geom_vline(data = sum2, aes(xintercept = mean, color = name), linetype = 2) +
  geom_text(data = sum2, aes(x = mean, y = length(unique(col_sum$sample)) + 0.75, label = paste(comma(mean / 1e6), "Mbp"))) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.75) +
  geom_text(aes(label = paste(comma(value / 1e6), "Mbp")), hjust = 0, vjust = 1) +
  scale_x_continuous(labels = comma) +
  scale_fill_manual(values = c("blue", "orange")) +
  scale_color_manual(values = c("blue", "orange")) +
  theme_cowplot() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_y_discrete(guide = guide_axis(n.dodge = 1)) +
  coord_cartesian(clip = "off") +
  ggtitle("The number of collapsed and predicted unassembled bases in HPRC HiFi assemblies") +
  ylab("") +
  xlab("Number of collapsed bp")
col.plot
my_ggsave("{odir}/collapses.pdf", height = 8, width = 12, plot = col.plot)
