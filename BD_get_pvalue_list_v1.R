#!/usr/bin/env Rscript

args <- commandArgs({trailingOnly=TRUE})
if (length(args)==0){
  stop("Zero arguments provided (two expected).n", call.=FALSE)
}

library("DescTools")

df <- read.csv(args[1], header=TRUE)
num_windows <- nrow(df)
pvalues <- vector(mode="numeric", length=num_windows) #Initialize vector (mem saver)

#### Calculate BD statistic and pvalue for each window comparison
for (i in 1:num_windows){
  ## Cell values with Fisher adjustment to prevent divide by zero occurrences
  bd_1 <- df$bd1[i] + 0.5
  bd_2 <- df$bd2[i] + 0.5
  bd_3 <- df$bd3[i] + 0.5
  bd_4 <- df$bd4[i] + 0.5
  bd_5 <- df$bd5[i] + 0.5
  bd_6 <- df$bd6[i] + 0.5
  bd_7 <- df$bd7[i] + 0.5
  bd_8 <- df$bd8[i] + 0.5
  ## Calculate Odds 
  odds_one <- ((bd_1) / (bd_1 + bd_2))
  odds_two <- ((bd_3) / (bd_3 + bd_4))
  odds_three <- ((bd_5) / (bd_5 + bd_6))
  odds_four <- ((bd_7) / (bd_7 + bd_8))
  max_odds <- max(c(odds_one, odds_two, odds_three, odds_four))
  if (max_odds %in% c(odds_two, odds_three, odds_four) == TRUE){
    #Skip calculation if max odds are not odds 1
    pvalues[i] <- 999
  } else {
    ## Calculate Breslow Day test
    bd_table <- xtabs(freq ~ ., cbind(expand.grid(phenotype=c("sterile", "fertile"),
                                                  window_one=c("focal", "non-focal"),
                                                  window_two=c("focal", "non-focal")),
                                      freq=c(bd_1, bd_2, bd_3, bd_4, bd_5, bd_6, bd_7, bd_8)))
    y <- BreslowDayTest(bd_table)$p.value
    pvalues[i] <- y
  }
  #sanity removes
  rm(bd_1, bd_2, bd_3, bd_4, bd_5, bd_6, bd_7, bd_8)
  rm(odds_one, odds_two, odds_three, odds_four, max_odds)
}

## Sanity matching test check
if (any(is.logical(pvalues) == TRUE)){
  writeLines('     - Warning: The number of tests performed did not match the number of window comparisons')
}

###Print to output
out_name <- paste(as.character(args[1], '_pvalue_list.csv',sep=''))
write.csv(pvalues, out_name, row.names=FALSE)
