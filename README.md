This repository contains the scripts used for the following study:

## Combinatorial expression of neurexin genes regulates glomerular targeting by olfactory sensory neurons

### Sung Jin Park<sup>1</sup>, I-Hao Wang<sup>1</sup>, Namgyu Lee<sup>2</sup>, Hao-Ching Jiang<sup>1</sup>, Takeshi Uemura<sup>3,4</sup>, Kensuke Futai<sup>5</sup>, Dohoon Kim<sup>6</sup>, Evan Z. Macosko<sup>7,8</sup>, Paul L. Greer<sup>1,*</sup>

<sup>1</sup>Program in Molecular Medicine, University of Massachusetts Medical School, Worcester, MA, USA,
<sup>2</sup>Department of Biomedical Science & Engineering, Dankook University, Cheonan, South Korea
<sup>3</sup>Division of Gene Research, Research Center for Advanced Science, Shinshu University, Nagano, Japan
<sup>4</sup>Institute for Biomedical Sciences, Interdisciplinary Cluster for Cutting Edge Research, Shinshu University, Nagano, Japan
<sup>5</sup>Department of Neurobiology and Brudnick Neuropsychiatric Research Institute, University of Massachusetts Medical School, Worcester, MA, USA
<sup>6</sup>Department of Molecular, Cell and Cancer Biology, University of Massachusetts Chan Medical School, Worcester, MA, USA
<sup>7</sup>Broad Institute of Harvard and MIT, Cambridge, MA, USA
<sup>9</sup>Department of Psychiatry, Massachusetts General Hospital, Boston, MA, USA, 
<sup>10</sup>These authors contributed equally, 
<sup>*</sup>Corresponding author.

```R
sessionInfo()
R version 4.3.0 (2023-04-21)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS 14.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggrepel_0.9.4        dplyr_1.1.4          GO.db_3.17.0         org.Mm.eg.db_3.17.0 
 [5] AnnotationDbi_1.62.2 IRanges_2.34.1       S4Vectors_0.38.2     Biobase_2.60.0      
 [9] BiocGenerics_0.46.0  reshape2_1.4.4       pheatmap_1.0.12      viridis_0.6.4       
[13] viridisLite_0.4.2    ggplot2_3.4.4        patchwork_1.1.3      harmony_1.2.0       
[17] Rcpp_1.0.11          Seurat_5.0.1         SeuratObject_5.0.1   sp_2.1-2            

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3      rstudioapi_0.15.0       jsonlite_1.8.8         
  [4] magrittr_2.0.3          ggbeeswarm_0.7.2        spatstat.utils_3.0-4   
  [7] farver_2.1.1            zlibbioc_1.46.0         vctrs_0.6.5            
 [10] ROCR_1.0-11             memoise_2.0.1           spatstat.explore_3.2-5 
 [13] RCurl_1.98-1.13         htmltools_0.5.7         sctransform_0.4.1      
 [16] parallelly_1.36.0       KernSmooth_2.23-22      htmlwidgets_1.6.4      
 [19] ica_1.0-3               plyr_1.8.9              plotly_4.10.3          
 [22] zoo_1.8-12              cachem_1.0.8            igraph_1.6.0           
 [25] mime_0.12               lifecycle_1.0.4         pkgconfig_2.0.3        
 [28] Matrix_1.6-4            R6_2.5.1                fastmap_1.1.1          
 [31] GenomeInfoDbData_1.2.10 fitdistrplus_1.1-11     future_1.33.0          
 [34] shiny_1.8.0             digest_0.6.33           colorspace_2.1-0       
 [37] tensor_1.5              RSpectra_0.16-1         irlba_2.3.5.1          
 [40] RSQLite_2.3.4           labeling_0.4.3          progressr_0.14.0       
 [43] fansi_1.0.6             spatstat.sparse_3.0-3   httr_1.4.7             
 [46] polyclip_1.10-6         abind_1.4-5             compiler_4.3.0         
 [49] withr_2.5.2             bit64_4.0.5             DBI_1.1.3              
 [52] fastDummies_1.7.3       MASS_7.3-60             tools_4.3.0            
 [55] vipor_0.4.5             lmtest_0.9-40           beeswarm_0.4.0         
 [58] httpuv_1.6.13           future.apply_1.11.0     goftest_1.2-3          
 [61] glue_1.6.2              nlme_3.1-164            promises_1.2.1         
 [64] grid_4.3.0              Rtsne_0.17              cluster_2.1.6          
 [67] generics_0.1.3          gtable_0.3.4            spatstat.data_3.0-3    
 [70] tidyr_1.3.0             data.table_1.14.10      utf8_1.2.4             
 [73] XVector_0.40.0          spatstat.geom_3.2-7     RcppAnnoy_0.0.21       
 [76] RANN_2.6.1              pillar_1.9.0            stringr_1.5.1          
 [79] spam_2.10-0             RcppHNSW_0.5.0          later_1.3.2            
 [82] splines_4.3.0           lattice_0.22-5          survival_3.5-7         
 [85] bit_4.0.5               deldir_2.0-2            tidyselect_1.2.0       
 [88] Biostrings_2.68.1       miniUI_0.1.1.1          pbapply_1.7-2          
 [91] gridExtra_2.3           scattermore_1.2         matrixStats_1.2.0      
 [94] stringi_1.8.3           lazyeval_0.2.2          codetools_0.2-19       
 [97] tibble_3.2.1            cli_3.6.2               uwot_0.1.16            
[100] xtable_1.8-4            reticulate_1.34.0       munsell_0.5.0          
[103] GenomeInfoDb_1.36.4     globals_0.16.2          spatstat.random_3.2-2  
[106] png_0.1-8               ggrastr_1.0.2           parallel_4.3.0         
[109] ellipsis_0.3.2          blob_1.2.4              dotCall64_1.1-1        
[112] bitops_1.0-7            listenv_0.9.0           scales_1.3.0           
[115] ggridges_0.5.5          leiden_0.4.3.1          purrr_1.0.2            
[118] crayon_1.5.2            rlang_1.1.2             cowplot_1.1.2          
[121] KEGGREST_1.40.1  
```
