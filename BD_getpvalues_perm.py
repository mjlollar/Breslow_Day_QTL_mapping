#!/usr/bin/env Rscript

### Author: Matt
### mjlollar1@gmail.com
### Last update: Aug 3rd 2023

## Input file, BD cell counts file
## one arguments required, input file)
## e.g. $ Rscript myinput_bdgroups_perm.csv


args <- commandArgs({trailingOnly=TRUE})
if (length(args)==0){
  stop("Zero arguments provided (two expected).n", call.=FALSE)
}

library("DescTools") #Library containing Breslow-Day Test

df <- read.csv(args[1], sep=',') #Expects comma-delim input file
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
    #print(y)
  }
  #sanity removes
  rm(bd_1, bd_2, bd_3, bd_4, bd_5, bd_6, bd_7, bd_8)
  rm(odds_one, odds_two, odds_three, odds_four, max_odds, bd_table, i, y)
}

## Sanity matching test check
if (any(is.logical(pvalues) == TRUE)){
  writeLines('     - Warning: The number of tests performed did not match the number of window comparisons')
}

# Divide pvalues and get lowest X-A and A-A null pvalue
pvalues_x <- pvalues[1:1108500] #Set of X-2 and X-3 comparisons
pvalues_a <- pvalues[1108500:num_windows] #2-3 comparisons
lowest_pvalue_x <- as.character(min(pvalues_x))
lowest_pvalue_a <- as.character(min(pvalues_a))

# Output lowest X-A and lowest A-A pvalue as separate txts
out_name_x <- paste(as.character(args[1]), '_x_null_pvalue.txt',sep='')
out_file_x <- file(out_name_x)
writeLines(lowest_pvalue, out_file_x)
out_name_a <- paste(as.character(args[1]), '_a_null_pvalue.txt',sep='')
out_file_a <- file(out_name_a)
writeLines(lowest_pvalue_a, out_file_a)

close(out_file_x)
close(out_file_a)
