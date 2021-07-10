if(! require("tidyverse")) install.packages("tidyverse")
if(! require("ggnewscale")) install.packages("ggnewscale")
#if(! require("ggrepel")) remove.packages("ggrepel"); install.packages("ggrepel"); #install.packages("https://cran.r-project.org/bin/macosx/contrib/4.0/ggrepel_0.9.1.tgz", repos = NULL, type="source")
library(ggrepel)
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
library(ggforce)
library(sys)
#install.packages("ggdist")
library(ggdist)
sampleinfo = fread("~/Desktop/EichlerVolumes/chm13_t2t/nobackups/assemblies_for_anlysis/sample_info/Master.tbl")
tmp1 = fread("~/Desktop/sampled-full_pangenie-38_genotyping-lpa_all.tsv"); tag="all"
tmp2 = fread("~/Desktop/sampled-full_pangenie-38_genotyping-variable20_all.tsv"); tag="20_invariable"
tmp3 = fread("~/Desktop/sampled-full_pangenie-38_genotyping-variable90_all.tsv"); tag="90_invariable"

# colors to use 
GRAY = "#2F4F4F"	
RED = "#af0404"
BLUE = "#3282b8"
NEWCOLOR = RED
OLDCOLOR = GRAY 

my_fread = function(file){
  df=fread(file)
  df$file=gsub("merged_sampled-full_pangenie-38_genotyping-|_all.tsv", "", basename(file))
  df %>%
    separate(file, into=c("Region", "Filter"), sep="_") %>%
    mutate(Filter=factor(recode(Filter,
                      `conf0`="None",
                      `conf1`="Lenient",
                      `conf4`="Strict"
                      ), levels = c("None","Lenient","Strict"))
           ) %>%
    merge(unique(sampleinfo[,c("sample","SuperPop")]), by=c("sample")) %>%
    data.table()
}

tmp = bind_rows(lapply(Sys.glob("~/Desktop/LPA_genotyping_tables/*.tsv"), my_fread))
tmp$Region = gsub("variable", "non-variable, n=", tmp$Region)
#tmp = bind_rows(list(all=tmp1, `20`=tmp2, `90`=tmp3), .id="tag")
#tmp$tag = factor(tmp$tag, levels = unique(tmp$tag))
df = tmp %>% 
  arrange(SuperPop) %>%
  mutate(sample = factor(sample, levels=unique(sample))) %>%
  data.table()
head(df)
dim(df)
colnames(df)
unique(df$Region)

for(z in unique(df$Region)){
p = ggplot(data=df[df$Region==z], aes(x=sample, y=genotype_concordance, group=Region)) +
  ggrepel::geom_text_repel(
    aes(
      label=paste(`nr_correct`, `nr_typed`, sep="/")
      ),
    direction = "y", nudge_y = 10
    ) + 
  ggrepel::geom_text_repel(
    aes(y=`non-ref_concordance`, 
        label=paste(`nr_correct_non-ref`, `nr_typed_non-ref`, sep="/")
        ),
    alpha=0.5,
    point.size = 5,
    direction = "y", nudge_y = -3
    ) + 
  geom_point(aes(color=SuperPop), size=3) +
  geom_point(aes(y=`non-ref_concordance`,color=SuperPop), size=3) +
  geom_segment(aes(xend=sample, yend=`non-ref_concordance`), size = 1, linetype=2, alpha=0.5, color="gray") +
  theme_minimal_hgrid() +
  scale_x_discrete(guide = guide_axis(n.dodge = 3)) + scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme(legend.position = "top") +
  facet_grid(Filter~.) +
  ylab("PanGenie genotyping concordance in LPA") +
  scale_color_brewer(palette="Set1");p
scale=1.5
ggsave(glue("~/Desktop/LPA_genotyping_tables/LPA_genotype_{z}.pdf"), height = 8*scale, width=12*scale, plot=p)
print(z)
}

df.frac = df %>% 
  mutate(  
    ref = (nr_correct - `nr_correct_non-ref`)/(nr_typed-`nr_typed_non-ref`)*100,
    non_ref = (`nr_correct_non-ref`)/(`nr_typed_non-ref`)*100
    ) 
summ = df.frac %>%
  group_by(Region, Filter) %>%
  summarise(
    mean_all = mean(genotype_concordance),
    median_all= median(genotype_concordance),
    sd_all = sd(genotype_concordance),
    mean_ref = mean(ref),
    median_ref = median(ref),
    sd_ref = sd(ref),
    mean_nonref = mean(non_ref),
    median_nonref =median(non_ref), 
    sd_nonref=sd(non_ref),
    n=n()
  ) %>% data.table()

gg.summ = summ %>% 
  pivot_longer(colnames(summ)[3:(length(summ)-1)]) %>%
  filter(name!="n")%>%
  separate(name, into = c("Measure","Subset"),sep="_") %>%
  data.table()

gg.summ

ggplot(data=gg.summ, aes(y=value, x=Filter, color=Subset))+
  geom_point() +
  geom_text_repel(aes(label=round(value,2)),
                  nudge_x = 0.5,
                  direction = "y"
                  )+
  facet_grid(Measure~Region, scales="free_y")+
  theme_minimal_grid()+
  theme(legend.position = "top")

summ
summ.plot = df.frac %>% pivot_longer(cols=c("genotype_concordance", "ref","non_ref")) %>%
  mutate(name=gsub("genotype_concordance", "all", name)) %>%
  mutate(name=gsub("non_ref", "non-ref", name)) %>%
  mutate(name=factor(name, levels=c("all", "ref", "non-ref"))) %>%
  ggplot(aes(x=Region, y=value, color=Filter)) +
  facet_col(~name, scales = "free")+#, space = "free") +
  theme_cowplot()+
  theme(legend.position = "top") + 
  geom_boxplot(
    width = .2, 
    ## remove outliers
    outlier.color = NA ## `outlier.shape = NA` works as well
  ) +
  ## add dot plots from {ggdist} package
  ggdist::stat_dots(
    aes(fill=Filter), alpha=0.6,
    ## orientation to the left
    side = "left", 
    ## move geom to the left
    justification = 1.15,
    width=0.8
  )+
  scale_color_manual(values = c(RED, GRAY, BLUE))+
  scale_fill_manual(values = c(RED, GRAY, BLUE))+
  ylab("Genotyping accuracy")+xlab("Subsets"); summ.plot
ggsave("~/Desktop/LPA_genotyping_tables/LPA_summ.pdf", plot = summ.plot, height = 8, width = 12)
