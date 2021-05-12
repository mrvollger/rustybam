#!/usr/bin/env Rscript
options(repos=structure(c(CRAN="http://cran.us.r-project.org")))
.libPaths(c("~/local/R/library", .libPaths()))

if(! require("tidyverse")) install.packages("tidyverse")
if(! require("ggnewscale")) install.packages("ggnewscale")
if(! require("ggrepel")) install.packages("ggrepel")
if(! require("data.table")) install.packages("data.table")
if(! require("glue")) install.packages("glue")
if(! require("RColorBrewer")) install.packages("RColorBrewer")
if(! require("scales")) install.packages("scales")
if(! require("cowplot")) install.packages("cowplot")
if(! require("argparse")) install.packages("argparse")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(! require("karyoploteR")) BiocManager::install("karyoploteR")
if(! require("GenomicRanges")) BiocManager::install("GenomicRanges")



# create parser object
indir="~/Desktop/EichlerVolumes/chm13_t2t/nobackups/assembly_alignments/rustybam/reference_alignment/bed/"
parser <- ArgumentParser()
parser$add_argument("-a", "--asm",  help="bed file with all the asm mapping", default = glue("{indir}/HG00733_1.bed"))
parser$add_argument("-b", "--asm2",  help="bed file with a second asm mapping", default = glue("{indir}/HG00733_2.bed"))
parser$add_argument("-k", "--karyotype",  help="karyotpye file for different genomes")
parser$add_argument("-p", "--plot",  help="output plot, must have .pdf ext.", default = "~/Desktop/ideogram.pdf")
args <- parser$parse_args()


asmdf<- function(filename, colors){
  asmvshg = read.table(filename, header=F)
  names(asmvshg) = c("chr", "start", "end", "name", "xqual")
  curcolor = 1
  lencolors = length(colors)
  precontig = ""
  asmcolor = NULL
  y = NULL
  for(i in 1:nrow(asmvshg) ){
    contig = as.character(asmvshg$name[i])
    if(contig != precontig){
      curcolor = (curcolor + 1) %% lencolors 
      precontig = contig
    }
    asmcolor = c(asmcolor, colors[curcolor + 1])
    y = c(y, curcolor/4)
  }
  asmvshg$color = asmcolor
  asmvshg$y = y
  asmvshg$y1 = asmvshg$y + .25
  return(asmvshg)
}

asmvshg = asmdf(args$asm,  c("#2081f9", "#f99820") ) 

if(!is.null(args$asm2)){
  asmvshg2 = asmdf(args$asm2,  c("#159934", "#99157a") )
}


cex = 0.5 

print("Plotting") 

pdf(file=args$plot, width = 9, height =11 )

if(is.null(args$asm2)){
  kp <- plotKaryotype(genome=GENOME, cytobands = CYTOFULL, chromosomes = NOM)
} else {
kp <- plotKaryotype(genome=GENOME, cytobands = CYTOFULL, chromosomes = NOM, plot.type = 2)
}

# adding asm bed number one
kpRect(kp, chr=asmvshg$chr, x0=asmvshg$start, x1=asmvshg$end, y0=asmvshg$y, y1=asmvshg$y1, col=asmvshg$color)

# adding second asm if there
if(!is.null(args$asm2)){
  kpRect(kp, chr=asmvshg2$chr, x0=asmvshg2$start, x1=asmvshg2$end, y0=asmvshg2$y, y1=asmvshg2$y1, col=asmvshg2$color, data.panel = 2)
}

dev.off()

