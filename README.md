# evapR
An R package to calculate evaporation from water surface 

## Installation

Package can be installed by:

```
devtools::install_github("tgmwri/evapR")
```

## Example

``` r
require(evapR)
require(dplyr)

options(stringsAsFactors = FALSE) #strongly recommended

dta <- readRDS(file.path(.dir, "data" , "radiace.rds"))

tbl <- compute_coef_table(dta)

eq_list <- tbl$ult_tab %>% arrange(desc(KGE)) %>% top_n(30, KGE) 
eq_list <- eq_list$final_eq
all_coef <- tbl$all_coef


eq_select <- find_fun(c("R", "Ta"), eq_list)
result <- apply_fun(dta, all_coef, eq_select)
```