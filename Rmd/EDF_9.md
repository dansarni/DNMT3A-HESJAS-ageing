Sarni et al., Extended Data Figure 9
================
dsarni
03-03-2026

## Extended Figure 9. Pax5 DNA hypermethylation is selectively lost in may B lymphocytes that progress to later differentiation stages

1.  Libraries used in this figure.

``` r
library(ggplot2)
```

2.  Import data

``` r
# Get heatmap file names
files.heatmap <- list.files(
      #path = "U:/Datastore/IGMM/ajackson-wrkgrp/People/Dan/DNMT3A_Project/EM/EM_PCR/EM_conversion/Nature_Genetics_Figures_Feb2026/files_for_figures",
      path = "../data/Figure_5_EDF9",
      pattern = "^heatmap.*\\.tsv$",
      full.names = TRUE
)

# Import files
for (i in seq_along(files.heatmap)) {
  assign(
    tools::file_path_sans_ext(basename(files.heatmap[i])),
    read.delim(files.heatmap[i]),
    envir = .GlobalEnv
  )
}

### Import violin plot df ###
# Pax5
violin_pax5_a_df <- read.table("../data/Figure_5_EDF9/Pax5_A_violin_df_03-03-2026.tsv", header = T)
violin_pax5_b_df <- read.table("../data/Figure_5_EDF9/Pax5_B_violin_df_03-03-2026.tsv", header = T)
violin_pax5_cf_df <- read.table("../data/Figure_5_EDF9/Pax5_CF_violin_df_03-03-2026.tsv", header = T)

# Hoxc13
violin_hoxc13_a_df <- read.table("../data/Figure_5_EDF9/Hoxc13_A_violin_df_03-03-2026.tsv", header = T)
violin_hoxc13_b_df <- read.table("../data/Figure_5_EDF9/Hoxc13_B_violin_df_03-03-2026.tsv", header = T)
violin_hoxc13_cf_df <- read.table("../data/Figure_5_EDF9/Hoxc13_CF_violin_df_03-03-2026.tsv", header = T)

# Hoxd13
violin_hoxd13_a_df <- read.table("../data/Figure_5_EDF9/Hoxd13_A_violin_df_03-03-2026.tsv", header = T)
violin_hoxd13_b_df <- read.table("../data/Figure_5_EDF9/Hoxd13_B_violin_df_03-03-2026.tsv", header = T)
violin_hoxd13_cf_df <- read.table("../data/Figure_5_EDF9/Hoxd13_CF_violin_df_03-03-2026.tsv", header = T)
```

3.  List all heatmap files

``` r
list_of_heatmap_df <- ls(pattern = "^heatmap") # all

list_of_heatmap_pax5 <- ls(pattern = "^heatmap_Pax5") # Pax5
list_of_heatmap_pax5
```

    ##  [1] "heatmap_Pax5_4k946_A_df_03-03-2026" "heatmap_Pax5_4k947_A_df_03-03-2026"
    ##  [3] "heatmap_Pax5_4k947_Z_df_03-03-2026" "heatmap_Pax5_4k948_A_df_03-03-2026"
    ##  [5] "heatmap_Pax5_4k948_Z_df_03-03-2026" "heatmap_Pax5_4k949_A_df_03-03-2026"
    ##  [7] "heatmap_Pax5_4k949_Z_df_03-03-2026" "heatmap_Pax5_B_4k1056_03-03-2026"  
    ##  [9] "heatmap_Pax5_B_4k1057_03-03-2026"   "heatmap_Pax5_B_4k1064_03-03-2026"  
    ## [11] "heatmap_Pax5_B_4k1065_03-03-2026"   "heatmap_Pax5_CF_4k1056_03-03-2026" 
    ## [13] "heatmap_Pax5_CF_4k1057_03-03-2026"  "heatmap_Pax5_CF_4k1064_03-03-2026" 
    ## [15] "heatmap_Pax5_CF_4k1065_03-03-2026"

``` r
list_of_heatmap_hoxc13 <- ls(pattern = "^heatmap_Hoxc13") # Hoxc13
list_of_heatmap_hoxc13
```

    ##  [1] "heatmap_Hoxc13_A_4k946_03-03-2026"   "heatmap_Hoxc13_A_4k947_03-03-2026"  
    ##  [3] "heatmap_Hoxc13_A_4k948_03-03-2026"   "heatmap_Hoxc13_A_4k949_03-03-2026"  
    ##  [5] "heatmap_Hoxc13_B_4k1056_03-03-2026"  "heatmap_Hoxc13_B_4k1057_03-03-2026" 
    ##  [7] "heatmap_Hoxc13_B_4k1064_03-03-2026"  "heatmap_Hoxc13_B_4k1065_03-03-2026" 
    ##  [9] "heatmap_Hoxc13_CF_4k1056_03-03-2026" "heatmap_Hoxc13_CF_4k1057_03-03-2026"
    ## [11] "heatmap_Hoxc13_CF_4k1064_03-03-2026" "heatmap_Hoxc13_CF_4k1065_03-03-2026"

``` r
list_of_heatmap_hoxd13 <- ls(pattern = "^heatmap_Hoxd13") # Hoxd13
list_of_heatmap_hoxd13
```

    ##  [1] "heatmap_Hoxd13_A_4k946_03-03-2026"   "heatmap_Hoxd13_A_4k947_03-03-2026"  
    ##  [3] "heatmap_Hoxd13_A_4k948_03-03-2026"   "heatmap_Hoxd13_A_4k949_03-03-2026"  
    ##  [5] "heatmap_Hoxd13_B_4k1056_03-03-2026"  "heatmap_Hoxd13_B_4k1057_03-03-2026" 
    ##  [7] "heatmap_Hoxd13_B_4k1064_03-03-2026"  "heatmap_Hoxd13_B_4k1065_03-03-2026" 
    ##  [9] "heatmap_Hoxd13_CF_4k1056_03-03-2026" "heatmap_Hoxd13_CF_4k1057_03-03-2026"
    ## [11] "heatmap_Hoxd13_CF_4k1064_03-03-2026" "heatmap_Hoxd13_CF_4k1065_03-03-2026"

4.  A function to plot heatmaps

``` r
plot_heatmap <- function(df, n, title = deparse(substitute(df))) {
  
  heatmap(
    as.matrix(df[, paste0("pos_", 1:n)]),
    Rowv = NA,
    Colv = NA,
    col = c("white", "lightgrey", "darkblue"),
    scale = "none",
    labRow = NA,
    main = title,
    cex.main = 0.7,
    labCol = 1:n
  )
}
```

5.  Plot heatmaps

- Pax5

``` r
for (name in list_of_heatmap_pax5) {
  plot_heatmap(get(name), 34, title = name)
}
```

![](EDF_9_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-6.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-7.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-8.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-9.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-10.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-11.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-12.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-13.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-14.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-5-15.png)<!-- -->

- Hoxc13

``` r
for (name in list_of_heatmap_hoxc13) {
  plot_heatmap(get(name), 70, title = name)
}
```

![](EDF_9_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-6-5.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-6-6.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-6-7.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-6-8.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-6-9.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-6-10.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-6-11.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-6-12.png)<!-- -->

- Hoxd13

``` r
for (name in list_of_heatmap_hoxd13) {
  plot_heatmap(get(name), 26, title = name)
}
```

![](EDF_9_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-7-5.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-7-6.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-7-7.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-7-8.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-7-9.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-7-10.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-7-11.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-7-12.png)<!-- -->

6.  Violin plots

List violin objects

``` r
list_of_violin <- ls(pattern = "^violin") # Pax5
list_of_violin
```

    ## [1] "violin_hoxc13_a_df"  "violin_hoxc13_b_df"  "violin_hoxc13_cf_df"
    ## [4] "violin_hoxd13_a_df"  "violin_hoxd13_b_df"  "violin_hoxd13_cf_df"
    ## [7] "violin_pax5_a_df"    "violin_pax5_b_df"    "violin_pax5_cf_df"

7.  A function for violin plots

``` r
plot_violin <- function(df, title = deparse(substitute(df))) {
  
v_plot <- ggplot(df, aes(x = animal, y = perc.meth, fill = geno))+
                geom_violin(scale = "width", adjust = 1.5)+
                stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", fatten = 1.5)+
                scale_fill_manual(values = c("+/+" = "#8c97ce", "W326R/+" = "#ec7c7e"))+
                labs(x = "", y = "% of DNA methylation", title = name) +
                ylim(c(0,100))+
                theme_classic() +
                theme(
                  strip.text = element_text(size = 12),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank()  
                )
print(v_plot)
}
```

8.  Plot violin plots

``` r
for (name in list_of_violin) {
  plot_violin(get(name), title = name)
}
```

![](EDF_9_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-10-5.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-10-6.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-10-7.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-10-8.png)<!-- -->![](EDF_9_files/figure-gfm/unnamed-chunk-10-9.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 4.5.0 (2025-04-11 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26100)
    ## 
    ## Matrix products: default
    ##   LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.utf8 
    ## [2] LC_CTYPE=English_United Kingdom.utf8   
    ## [3] LC_MONETARY=English_United Kingdom.utf8
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.utf8    
    ## 
    ## time zone: Europe/London
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_3.5.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.6.5        cli_3.6.5          knitr_1.50         rlang_1.1.6       
    ##  [5] xfun_0.52          generics_0.1.4     labeling_0.4.3     glue_1.8.0        
    ##  [9] htmltools_0.5.8.1  scales_1.4.0       rmarkdown_2.29     grid_4.5.0        
    ## [13] evaluate_1.0.4     tibble_3.3.0       fastmap_1.2.0      yaml_2.3.10       
    ## [17] lifecycle_1.0.4    compiler_4.5.0     dplyr_1.1.4        RColorBrewer_1.1-3
    ## [21] pkgconfig_2.0.3    rstudioapi_0.17.1  farver_2.1.2       digest_0.6.37     
    ## [25] R6_2.6.1           tidyselect_1.2.1   pillar_1.11.0      magrittr_2.0.3    
    ## [29] withr_3.0.2        tools_4.5.0        gtable_0.3.6
