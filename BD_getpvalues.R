#!/usr/bin/env Rscript

###Author: Matt Lollar
### mjlollar1@gmail.com
### Last update: Aug 3rd 2023

## Input file, BD cell counts file
## Two arguments are required, input file and type of scan (forward (0) or reverse (2))
## e.g. $ Rscript myinput_bdgroups_reverse.csv 2

args <- commandArgs({trailingOnly=TRUE})
if (length(args)==0){
  stop("Zero arguments provided (two expected).n", call.=FALSE)
}

library("DescTools") #Library containing Breslow-Day Test

df <- read.csv("54_403_bdgroups_080223_forward_scan_bd_cells.csv", sep=',')
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

### Get windows compared
### Chromosome cutoffs are hard coded, fits 50kb F2 files (2579 total windows)
### Order of comparisons matters, in that p-values ordered as X-2,X-3, etc. for correct index match to p-value
if (num_windows < 3000) { #Catches unidirectional scans
  window_one <- seq(1,2579, by=1) # Uni comparison
  out_uni_df <- data.frame(window_one, pvalues)
  write.csv(out_uni_df, "pvalue_list_unidirectional_output.csv", row.names=FALSE)
} else{
  window_one <- vector(mode="numeric", length=num_windows)
  window_two <- vector(mode="numeric", length=num_windows)
  counter <- 1
  if (args[2] == 0){
    # X - 2
    for(win_1 in 1:545){
      for (win_2 in 546:1524){
        window_one[counter] <- win_1
        window_two[counter] <- win_2
        counter <- counter + 1
      }
    }
    rm(win_1, win_2) #sanity
    # X - 3
    for (win_1 in 1:545){
      for (win_2 in 1525:2579){
        window_one[counter] <- win_1
        window_two[counter] <- win_2
        counter <- counter + 1
      }
    }
    rm(win_1, win_2) #sanity
    # 2 - 3
    for (win_1 in 546:1524){
      for (win_2 in 1525:2579){
        window_one[counter] <- win_1
        window_two[counter] <- win_2
        counter <- counter + 1
      }
    }
    rm(win_1, win_2) #sanity
  } else {
    # 2 - X
    for(win_1 in 546:1524){
      for (win_2 in 1:545){
        window_one[counter] <- win_1
        window_two[counter] <- win_2
        counter <- counter + 1
      }
    }
    rm(win_1, win_2) #sanity
    # 3 - X
    for (win_1 in 1525:2579){
      for (win_2 in 1:545){
        window_one[counter] <- win_1
        window_two[counter] <- win_2
        counter <- counter + 1
      }
    }
    rm(win_1, win_2) #sanity
    # 3 - 2
    for (win_1 in 1525:2579){
      for (win_2 in 546:1524){
        window_one[counter] <- win_1
        window_two[counter] <- win_2
        counter <- counter + 1
      }
    }
  }
  out_bi_df <- data.frame(window_one, window_two, pvalues)
  write.csv(out_bi_df, "pvalue_list_bidirectional_output.csv", row.names=FALSE)
}
