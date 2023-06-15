Density plots of clustering performance and BIRCH radius on 4 simulated
SM datasets
================
2022-11-12

Load libraries

``` r
library(knitr)
library(ggplot2)
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.5.0 
    ## ✔ readr   2.1.3      ✔ forcats 0.5.2 
    ## ✔ purrr   1.0.1      
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(ggridges)
library(ggthemes)
library(tidyr)
library(reshape2)
```

    ## 
    ## Attaching package: 'reshape2'
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

``` r
library(ggbeeswarm)
library(ggpubr)
library(eoffice)
```

Process clustering results into evaluation metrics

``` r
setwd<-("~/Desktop/AnchorClustering_Code/cha2_birch_radius_4datasets")
INPUT_FILES <- list.files(path="~/Desktop/AnchorClustering_Code/cha2_birch_radius_4datasets", pattern="HDC_*", full.names=T,recursive = FALSE)

OUT_FILE <- 'parameter_perf_4samples_birch.sim.tsv'

for(file in INPUT_FILES) {
        
dbanc<-read.table(file,header=TRUE,sep="\t")

 db<-dbanc
  
  tpos <- db %>%
    group_by(truth, cluster, junction_length) %>%
    summarise(tpos=0.5*n()*(n()-1), tpos_nseq=n())
  
  # Sum true positive pairs by condition positive (actual clonal relatives)
  tpos_tot <- tpos %>%
    group_by(truth, junction_length) %>%
    summarise(tpos=sum(tpos))
  
  # Count condition positive and add column to tpos_tot
  cond_pos <- db %>%
    group_by(truth, junction_length) %>%
    summarise(cond_pos=0.5*n()*(n()-1))
  tpos_tot <- left_join(tpos_tot, cond_pos)
  
  # Summarise tpr (sensitivity) by junction length
  tpos_sum <- tpos_tot %>%
    group_by(junction_length) %>%
    summarise(tpos=sum(tpos), cond_pos=sum(cond_pos))
  
  # Count total possible pairwise relations within junction length to get cond_neg
  cond_tot <- db %>%
    group_by(junction_length) %>%
    summarise(cond_all=0.5*n()*(n()-1))
  tpos_sum <- left_join(tpos_sum, cond_tot)
  tpos_sum$cond_neg <- tpos_sum$cond_all - tpos_sum$cond_pos
  
  # Sum true positive pairs by test positive (inferred clonal relatives)
  ppv_tot <- tpos %>%
    group_by(cluster, junction_length) %>%
    summarise(tpos=sum(tpos))
  
  # Count test positive and add column to ppv_tot
  test_pos <- db %>%
    group_by(cluster, junction_length) %>%
    summarise(test_pos=0.5*n()*(n()-1))
  ppv_tot <- left_join(ppv_tot, test_pos)
  
  # Summarize ppv by junction length
  ppv_sum <- ppv_tot %>%
    group_by(junction_length) %>%
    summarise(tpos=sum(tpos), test_pos=sum(test_pos))
  
  tot_sum <- left_join(tpos_sum, ppv_sum)
  # False Positives = Test positive - True positives
  tot_sum$fpos <- tot_sum$test_pos - tot_sum$tpos
  # Rename junction_length to junction_length
  names(tot_sum) <- tolower(names(tot_sum))
  
  file<-gsub("/Users/haiyangchang/Desktop/AnchorClustering_Code/cha2_birch_radius_4datasets/","",file)
  dataset<-gsub("MS-","", file)
  dataset<-gsub("HDC_","", dataset)
  dataset2<-gsub("_sim-pass.tab.csv","", dataset)
  tot_sum$norm_mini<-unlist(strsplit(dataset2, split="_"))[3]
  tot_sum$anc_num<-unlist(strsplit(dataset2, split="_"))[2]
  tot_sum$cdr3_length<-unlist(strsplit(dataset2, split="_"))[1]
  tot_sum$radius<-unlist(strsplit(dataset2, split="_"))[4]
  tot_sum$sample<-paste0(unlist(strsplit(dataset2, split="_"))[5],sep = '-',unlist(strsplit(dataset2, split="_"))[6])

  
  # Write summary false negative table
  WRITE_COLS <- !file.exists(OUT_FILE)
  write.table(tot_sum, file=OUT_FILE, append=T, 
              quote=F, sep='\t', row.names=F, col.names=WRITE_COLS)
}
```

Summarize clustering performance metrics

``` r
q<-read.table("~/Desktop/AnchorClustering_Code/cha2_birch_radius_4datasets/parameter_perf_4samples_birch.sim.tsv",header = TRUE,sep ="\t")
q$sen<-100*(as.numeric(q$tpos)/as.numeric(q$cond_pos))
q$ppv<-100*(as.numeric(q$tpos)/as.numeric(q$test_pos))
q$spec<-100*(as.numeric(q$fpos)/as.numeric(q$cond_neg))
q$fm<-(2*as.numeric(q$sen)*as.numeric(q$ppv))/(as.numeric(q$sen)+as.numeric(q$ppv))
q<-q[complete.cases(q), ] 

q$data<-gsub("-.*","",q$sample)

q2<-q %>% distinct(junction_length, norm_mini, fm, anc_num, radius, sample)
```

Explore the correlation between BIRCH radius and minimum distance
ratio(showing anchor numbers of that junction length) on the clustering
performance

``` r
q2$radius<-as.factor(q2$radius)
`Minimum distance ratio`<-as.factor(q2$norm_mini)
q2$`Anchor number`<-q2$anc_num


ggplot(q2, aes(x = fm, y = radius, fill = `Minimum distance ratio`, point_size = `Anchor number`)) +
  geom_density_ridges(
    aes(point_shape = `Minimum distance ratio`), 
    alpha = .2, point_alpha = 2, jittered_points = TRUE) + scale_point_size_continuous(range  = c(0.5, 5), 
                        limits = c(2, 160), 
                        breaks = c(2, 10, 25, 50, 100))+  
  scale_point_color_hue(l = 40) +  facet_wrap(~sample)+
  scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23, 24, 25))+
  labs(x = "F-measure", 
       y = "BIRCH radius")+theme_few()+
                theme(axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")) 
```

    ## Picking joint bandwidth of 1.01

    ## Picking joint bandwidth of 1.57

    ## Picking joint bandwidth of 2.09

    ## Picking joint bandwidth of 1.27

![](radius_series_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Explore the correlation between BIRCH radius and minimum distance ratio
on the clustering performance

``` r
q2$radius<-as.factor(q2$radius)
`Minimum distance ratio`<-as.factor(q2$norm_mini)
q2$`Anchor number`<-q2$anc_num

ggplot(q2, aes(x = fm, y = radius, fill = `Minimum distance ratio`)) +
  geom_density_ridges(
    aes(point_shape = `Minimum distance ratio`), 
    alpha = .2, point_alpha = 2, jittered_points = TRUE) +
  scale_point_color_hue(l = 40) +  facet_wrap(~sample)+
  scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23, 24, 25))+
  labs(x = "F-measure", 
       y = "BIRCH radius")+theme_few()+
        theme(axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")) 
```

    ## Picking joint bandwidth of 1.01

    ## Picking joint bandwidth of 1.57

    ## Picking joint bandwidth of 2.09

    ## Picking joint bandwidth of 1.27

![](radius_series_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Explore the correlation between BIRCH radius and minimum distance ratio
on the clustering performance

``` r
q2$radius<-as.factor(q2$radius)
q2$`Minimum distance ratio`<-as.factor(q2$norm_mini)
q2$`anchor number`<-q2$anc_num
ggplot(q2, aes(x = fm, y = radius, fill = `Minimum distance ratio`)) +
  geom_density_ridges(alpha = .6,) +
  facet_wrap(~sample)+
  scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23, 24, 25))+
  labs(x = "F-measure", 
       y = "BIRCH radius")+theme_few()+
           theme(axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")) 
```

    ## Picking joint bandwidth of 1.01

    ## Picking joint bandwidth of 1.57

    ## Picking joint bandwidth of 2.09

    ## Picking joint bandwidth of 1.27

![](radius_series_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
#ggsave("sup_pagea.jpg", plot = last_plot(), device = "jpeg", dpi = 600, width = 10, height = 7)

#dev.off()
```

Zoom in the BIRCH radius as 0.5 and show the correlation with minimum
distance on the clustering performance.

``` r
q3<-q[q$radius=="0.5",]

`Minimum distance ratio`<-as.factor(q3$norm_mini)
q3$`Anchor number`<-q3$anc_num

ggplot(q3, aes(x = fm, y = `Minimum distance ratio`, fill = `Minimum distance ratio`,point_size = `Anchor number`)) +
  geom_density_ridges(
    aes(point_shape = `Minimum distance ratio`), 
    alpha = .2, point_alpha = 2, jittered_points = TRUE) + scale_point_size_continuous(range  = c(0.5, 5), 
                        limits = c(2, 160), 
                        breaks = c(2, 10, 25, 50, 100))+  
  scale_point_color_hue(l = 40) +  facet_wrap(~sample)+ 
  scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23, 24, 25))+
  labs(x = "F-measure", 
       y = "Minimum distance ratio")+theme_few()+
   theme(axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
```

    ## Picking joint bandwidth of 0.727

    ## Picking joint bandwidth of 0.255

    ## Picking joint bandwidth of 0.979

    ## Picking joint bandwidth of 0.893

![](radius_series_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
#ggsave("sup_pageb.jpg", plot = last_plot(), device = "jpeg", dpi = 600, width = 10, height = 7)

#dev.off()
```
