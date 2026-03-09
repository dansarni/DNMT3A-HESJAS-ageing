Sarni et al., Figure 2
================
dsarni
24-02-2026

## Figure 2. Dnmt3a w326r/+ mice model HESJAS syndrome and exhibit multiple age-associated phenotypes

1.  Libraries used in this figure.

``` r
library(dplyr)
library(ggplot2)
library(RColorBrewer)
```

2.  Import data - DNA hypermethylated DMRs in bone marrow from +/+ and
    w326r/+ ages: 4-day, 23-day, 3-month and 1-year old mice.

``` r
# Import data tables

bm <- read.csv("../data/Figure_2/Figure_2h_bm_all_dmr_mean_age.csv")
```

3.  Convert the mean age data into a format suitable for ggplot

``` r
bm.df <- data.frame(Time = rep(c(rep("4d", length(bm[[1]])),
                             rep("23d", length(bm[[1]])),
                             rep("3m", length(bm[[1]])),
                             rep("1yr", length(bm[[1]]))),2),
                    Geno = c(rep("+/+",length(bm[[1]])*4),
                             rep("W326R/+",length(bm[[1]])*4)),
                    Value = c(bm$c.4d,
                              bm$c.23d,
                              bm$c.3m,
                              bm$c.1yr,
                              bm$m.4d,
                              bm$m.23d,
                              bm$m.3m,
                              bm$m.1yr))

bm.df$Group <- interaction(bm.df$Geno, bm.df$Time, sep = "_")
```

4.  Order samples for plotting.

``` r
group_levels <- levels(factor(bm.df$Group))
group_levels <- c(group_levels[7],
                  group_levels[3],
                  group_levels[5],
                  group_levels[1],
                  group_levels[8],
                  group_levels[4],
                  group_levels[6],
                  group_levels[2])
group_levels
```

    ## [1] "+/+_4d"      "+/+_23d"     "+/+_3m"      "+/+_1yr"     "W326R/+_4d" 
    ## [6] "W326R/+_23d" "W326R/+_3m"  "W326R/+_1yr"

``` r
bm.df$Group <- factor(bm.df$Group, levels = group_levels)
```

### Figure 2.h. Boxplot for mean DNA methylation per age per genotype

``` r
ggplot(bm.df, aes(x = Geno, y = Value, fill = Group)) +
  stat_boxplot(
    position = position_dodge(0.8),
    geom = "errorbar",
    width = 0.4
  ) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               colour = NA) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.75,
    colour = "black",
    position = position_dodge(0.8)
  ) +
  scale_fill_manual(values = c(
    "+/+_4d" = "#d4def1",
    "+/+_23d" = "#90bde5",
    "+/+_3m" = "#516fb5",
    "+/+_1yr" = "#292c6a",
    "W326R/+_4d" = "#f7b0b2",
    "W326R/+_23d" = "#f17275",
    "W326R/+_3m" = "#ec2f2d",
    "W326R/+_1yr" = "#9c1d22"
  )) +
  labs(x = "Genotype", y = "B-Value", fill = "Group") +
  theme_classic()
```

![](Figure_2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

5.  Compute the *P-values* between ages per genotype

- wild-type mice

``` r
# Kruskal-Wallis for +/+
bm.df.wt <- bm.df[bm.df$Geno == "+/+",]
  
kruskal.test(Value~Time, data = bm.df.wt)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Value by Time
    ## Kruskal-Wallis chi-squared = 136.24, df = 3, p-value < 2.2e-16

``` r
pairwise.wilcox.test(bm.df.wt$Value, bm.df.wt$Time, p.adjust.methods = "BH")
```

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
    ## 
    ## data:  bm.df.wt$Value and bm.df.wt$Time 
    ## 
    ##     1yr     23d     3m     
    ## 23d 2.9e-08 -       -      
    ## 3m  0.31    4.7e-06 -      
    ## 4d  < 2e-16 1.4e-06 < 2e-16
    ## 
    ## P value adjustment method: holm

- W326R/+ mice

``` r
# Kruskal-Wallis for W326R/+
bm.df.mut <- bm.df[bm.df$Geno == "W326R/+",]
  
kruskal.test(Value~Time, data = bm.df.mut)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Value by Time
    ## Kruskal-Wallis chi-squared = 495.47, df = 3, p-value < 2.2e-16

``` r
pairwise.wilcox.test(bm.df.mut$Value, bm.df.mut$Time, p.adjust.methods = "BH")
```

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
    ## 
    ## data:  bm.df.mut$Value and bm.df.mut$Time 
    ## 
    ##     1yr     23d     3m     
    ## 23d < 2e-16 -       -      
    ## 3m  1.6e-05 2.0e-06 -      
    ## 4d  < 2e-16 < 2e-16 < 2e-16
    ## 
    ## P value adjustment method: holm

Compute the *P-values* between genotypes (day 4)

``` r
# Two-sided, paired Wilcoxon rank sum test for +/+ versus W326R/+
day4 <- subset(bm.df, Time == "4d")
pairwise.wilcox.test(day4$Value, day4$Geno, p.adjust.methods = "BH")
```

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
    ## 
    ## data:  day4$Value and day4$Geno 
    ## 
    ##         +/+   
    ## W326R/+ <2e-16
    ## 
    ## P value adjustment method: holm

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
    ## [1] RColorBrewer_1.1-3 ggplot2_3.5.2      dplyr_1.1.4       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.6.5       cli_3.6.5         knitr_1.50        rlang_1.1.6      
    ##  [5] xfun_0.52         generics_0.1.4    labeling_0.4.3    glue_1.8.0       
    ##  [9] htmltools_0.5.8.1 scales_1.4.0      rmarkdown_2.29    grid_4.5.0       
    ## [13] evaluate_1.0.4    tibble_3.3.0      fastmap_1.2.0     yaml_2.3.10      
    ## [17] lifecycle_1.0.4   compiler_4.5.0    pkgconfig_2.0.3   rstudioapi_0.17.1
    ## [21] farver_2.1.2      digest_0.6.37     R6_2.6.1          tidyselect_1.2.1 
    ## [25] pillar_1.11.0     magrittr_2.0.3    withr_3.0.2       tools_4.5.0      
    ## [29] gtable_0.3.6
