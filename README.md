# TarCA.beta

The method, termed as targeting coalescent analysis (TCA), computes for all cells of a tissue the average coalescent rate at the monophyletic clades of the target tissue, the inverse of which then measures the progenitor number of the tissue. Any predefined population could be investigated with TCA, independent of pre-set markers.

<p align="center">
  <img src="https://github.com/shadowdeng1994/TarCA/blob/main/inst/CoalescentTheory.png" width=65% height=65%>
</p>

__* To achieve a higher computational efficiency, we rewrite the entire package in a more compact way. :smile:__  

## System requirement
* Dependent packages: dplyr, tidyr, tibble, ggplot2, ggtree
* Require R (>= 3.5.0).

## Install

```
install.packages('devtools')
devtools::install_github('shadowdeng1994/TarCA.beta')
```
Installation would finish in about one minute. 

## Quickstart
```
library("TarCA")
```
- The following files are needed for TarCA.
1. A tree file of class "phylo" with node labels.
> ((Cell_1,((Cell_2,Cell_3)Node_4,(Cell_4,Cell_5)Node_5)Node_3)Node_2,(((Cell_6,Cell_7)Node_8,(Cell_8,Cell_9)Node_9)Node_7,Cell_10)Node_6)Node_1;
2. A dataframe with columns *TipLabel* and *TipAnn*, representing tip labels on the tree file and corresponding cell annotations.

> | TipLabel | TipAnn |
> | --- | --- |
> | Cell_1 | O1 |
> | Cell_2 | O1 |
> | Cell_3 | O1 |
> | Cell_4 | O2 |
> | Cell_5 | O2 |
> | Cell_6 | O2 |
> | Cell_7 | O3 |
> | Cell_8 | O3 |
> | Cell_9 | O3 |
> | Cell_10 | O3 |

- Effective number of progenitor can be inferred with `Np_Estimator`.
- Modified algorithm for detection of lineage specific expression upregulation (LEU) can be called with `LEU_Estimator`.
- (optional) All intermediate data are stored in ExTree format (control with _ReturnExTree_, default FALSE).

## Estimate Np with exemplar dataset.
- Load exemplar dataset. 
```
load(system.file("Exemplar","Exemplar_TCA.RData",package = "TarCA.beta"))
tmp.tree <- ExemplarData_1$Tree
tmp.ann <- ExemplarData_1$Ann
```
- Inferring Np with `Np_Estimator`.
```
tmp.result <- Np_Estimator(
  Tree = tmp.tree,
  Ann = tmp.ann
)
```
> ===> Checking input files.  
> ===> Converting to ExTree.  
> ===> Adding AllDescendants.  
> ===> Adding MonoClades.  
> ===> Estimating Np.  


- Then return a dataframe containing the Np estimation.  
> | TipAnn | MonoInfo | Total | Np |
> | --- | --- | --- | --- |
> | O0 | 1 (1), 2 (2) | 5 | 5 |
> | O1 | 1 (6), 2 (11), 3 (1), 5 (1) | 36 | 26.2 |
> | O2 | 1 (35), 2 (17), 3 (4), 4 (2), 8 (1) | 97 | 67.5 |
> | O3 | 1 (66), 2 (38), 3 (11), 4 (4), 5 (2), 7 (1) | 208 | 158 |
> | O4 | 1 (50), 2 (24), 3 (6), 4 (3), 5 (1) | 133 | 125 |
> | O5 | 1 (71), 2 (38), 3 (13), 4 (5) | 206 | 197 |
> | O6 | 1 (32), 2 (23), 3 (9), 7 (1) | 112 | 87.5 |
> | O7 | 1 (50), 2 (37), 3 (10), 4 (3), 6 (1) | 172 | 147 |
> | O8 | 1 (5), 2 (3) | 11 | 18.3 |
> | O9 | 1 (12), 2 (1), 3 (2) | 20 | 27.1 |

This process is estimated to be completed in about __3__ seconds.

## Detect LEU with exemplar dataset.
- Load exemplar dataset.
```
load(system.file("Exemplar","Exemplar_LEU.RData",package = "TarCA.beta"))
tmp.tree <- ExemplarData_2$Tree
tmp.ann <- ExemplarData_2$Ann
```
- Inferring Np with `LEU_Estimator`.
```
tmp.result <- LEU_Estimator(
  Tree = tmp.tree,
  Ann = tmp.ann
)
```
> ===> Checking input files.  
> ===> Converting to ExTree.  
> ===> Adding AllDescendants.  
> ===> Adding ExpreBias.  
> ===> Adding FilterBiasParent.  
> ===> Estimating Np.  

- Then return a dataframe containing the Np estimation.  
> | TipAnn | MonoInfo | Total | Np |
> | --- | --- | --- | --- |
> | TRUE | 1 (23), 2 (5), 4 (2) | 41 | 48.2 |

This process is estimated to be completed in about __2__ seconds.

## Contributing
### Contributors
Shanjun Deng, shadowdeng1994@gmail.com.

## Citations
When using TarCA please cite:
