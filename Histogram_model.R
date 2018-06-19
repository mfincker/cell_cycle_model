#!/usr/bin/env Rscript

# Histogram_model.R
# Similar to Histogram_model.Rmd but makes it command line friendly,
# and doesn't produce a pdf.
# Only outputs the fits file and best fit plot

# Command line call:
# first argument should be the full path to the histo.tsv file
# ! Update path the Histogram_model.Rmd if running on a different computer.
  

library(tidyverse)
library(stringr)

## HISTOGRAM FILE

args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop(paste("Supply the full path to the histo file.\n",
             "Usage: Rscript --vanilla Histogram_model.R <full input file path>"),
       call. = FALSE)
}

histo_file <- args[1]

## CONSTANTS

chr_size <- 1.3

doubling_times <- #doubling times in days
  tibble(reactor = c("A", "Z", "M", "R"),
         doubling_time = c(12.1, 38.8, 48.5, 23.1))

# ## TEST
# 
# histo_file <- "/Users/maeva/Research/FACS/carmen-20180405/model/M_rep1_sample17_pico_rnase1_histo.tsv"
# fits_file <- str_replace(histo_file, "tsv", "fits")
# histo_test <- 
#   histo_file %>% 
#   read_tsv()
# 
# curr_reactor <- 
#   histo_file %>% 
#   str_match("\\/(.)_rep") %>% 
#   .[2]
# 
# curr_tau <- 
#   doubling_times %>% 
#   dplyr::filter(reactor == curr_reactor) %>% 
#   .$doubling_time


## FUNCTIONS 

### Get main peak channel

getC1 <- function(data_binned) {
  # Given a df of binned DNA data (channel ~ count), returns the C1 channel
  c1_pico <- 
    data_binned %>% 
    top_n(5, count) %>% 
    .$pico %>% 
    mean()
  
  c1_bin <- 
    data_binned %>% 
    mutate(dif_c1 = abs(pico - c1_pico)) %>% 
    top_n(n = -1,  wt = dif_c1) %>% 
    .$n - 1
  
  theo_c1_pico <- 
    data_binned %>% 
    filter(n == c1_bin) %>% 
    .$pico
  
  c2_bin <- 
    data_binned %>% 
    mutate(dif_c2 = abs(pico - 2 * theo_c1_pico)) %>% 
    top_n(n = -1, wt = dif_c2) %>% 
    .$n
  
  return(list(c1_bin, c2_bin))
}

#### TEST
# getC1(histo_test)
# histo_test %>% 
#   ggplot(aes(n, count)) +
#   geom_line() +
#   geom_vline(xintercept = getC1(histo_test) %>% unlist())


### DNA theoretical distribution

#### Theoretical count in a given channel

theoretical_count <- function(bin, start_bin, end_bin, B, C, tau, total_n, time_per_channel) {
  # Given a c1 channel, B, C, tau and total counts, returns the theoretical counts
  # in a specific channel.
  
  if (bin < start_bin | bin > end_bin) {
    count_theo <- 0
  } else if (bin == start_bin) {
    count_theo <- total_n * pop_percent(0, B + time_per_channel/2, tau)
  } else if (bin == end_bin) {
    count_theo <- total_n * pop_percent(B + C - time_per_channel/2, tau, tau)
  } else {
    count_theo <- total_n * pop_percent(B + time_per_channel/2 + time_per_channel * (bin - start_bin - 1),
                                        B + time_per_channel/2 + time_per_channel * (bin - start_bin), 
                                        tau) 
  }
  
  return(count_theo)
}

### TEST
# theoretical_count(11, 10, 12, 0.74, 0.12, 1, 100, 0.12/2)
# theoretical_count(10, 10, 12, 0.74, 0.12, 1, 100, 0.12/2)
# theoretical_count(12, 10, 12, 0.74, 0.12, 1, 100, 0.12/2)
# 
# theoretical_count(11, 10, 12, 0.74, 0.12, 1, 100, 0.12/2) +
# theoretical_count(10, 10, 12, 0.74, 0.12, 1, 100, 0.12/2) +
# theoretical_count(12, 10, 12, 0.74, 0.12, 1, 100, 0.12/2)
# 
# theoretical_count(5, 5, 8, 0.74, 0.12, 1, 100, 0.12/3)
# theoretical_count(6, 5, 8, 0.74, 0.12, 1, 100, 0.12/3)
# theoretical_count(7, 5, 8, 0.74, 0.12, 1, 100, 0.12/3)
# theoretical_count(8, 5, 8, 0.74, 0.12, 1, 100, 0.12/3)
# theoretical_count(4, 5, 8, 0.74, 0.12, 1, 100, 0.12/3)
# theoretical_count(9, 5, 8, 0.74, 0.12, 1, 100, 0.12/3)
# 
# theoretical_count(5, 5, 8, 0.74, 0.12, 1, 100, 0.12/3) +
# theoretical_count(6, 5, 8, 0.74, 0.12, 1, 100, 0.12/3) +
# theoretical_count(7, 5, 8, 0.74, 0.12, 1, 100, 0.12/3) +
# theoretical_count(8, 5, 8, 0.74, 0.12, 1, 100, 0.12/3)


#### Percentage of cells in a given time slice

pop_percent <- function(t1, t2, t_total) {
  2 * (2 ^ (-t1 / t_total) - 2 ^ (-t2 / t_total))
}

### TEST
# pop_percent(0, 0.74, 1)
# pop_percent(0.86, 1, 1)
# pop_percent(0, 1, 1)

### DNA histogram model for replication in current generation (slow growth)

DNAhisto_currentGeneration <- 
  function(histo, B, C, tau, std_dev, std_dev_var, chr_size, plot = FALSE, filename = NULL) {
    # Input is :
    #       - B: duration of B phase
    #       - C: duration of C phase
    #       - tau: doubling time
    #       - std_dev: std of C1 peak
    #       - std_dev_var: increase of std_dev per channel
    #       - chr_size: chromosome size
    # Optional input:
    #       - plot: should the function plot the model (default: FALSE)
    #       - filename: output file name (default: '')
    
    # if ( B + C > tau) {
    #   stop("DNA replication didn't start in current generation.")
    # }
    
    peaks <- getC1(histo)
    c1_peak <- peaks[[1]]
    c2_peak <- peaks[[2]]
    total_count <- histo$count %>% sum()

    start_bin <- c1_peak
    end_bin <- c2_peak
    time_per_channel <- C / (end_bin - start_bin)

    
    histo <-
      histo %>%
      mutate(count_theo = map_dbl(n, 
                                  ~theoretical_count(., start_bin, end_bin, B, C, tau, total_count, time_per_channel))
      )
    
    print(histo %>% summarise_all(sum))
    
    #std_dev with incremental change
    histo <-
      histo %>%
      mutate(std_dev = std_dev + std_dev_var * min(max(0, row_number() - c1_peak), c2_peak),
             count_model = map2(n, std_dev,
                                function(x, y) {histo$count_theo * dnorm(histo$n, x, y) %>%
                                    as_tibble}) %>%
               bind_cols() %>%
               summarise_all(sum) %>%
               gather(key = x, value = y) %>%
               .$y)
    
    # # std_dev as 5% CV of the mean
    # histo <-
    #   histo %>%
    #   mutate(std_dev = 0.20 * n,
    #          count_model = map2(n, std_dev,
    #                             function(x, y) {histo$count_theo * dnorm(histo$n, x, y) %>%
    #                                 as_tibble}) %>%
    #            bind_cols() %>%
    #            summarise_all(sum) %>%
    #            gather(key = x, value = y) %>%
    #            .$y)

    deviation <-
      histo %>%
      mutate(dif = (sqrt(count) - sqrt(count_model))^2) %>%
      # filter(n > 0.5 * c1_peak) %>% # ignore the debris / cells without full DNA
      filter(n < 1.25 * c2_peak, n > 0.25 * c1_peak) %>%
      transmute(dif = dif) %>%
      sum() %>%
      # sqrt(.) / 255
      sqrt(.) / round(1.25 * c2_peak - 0.25 * c1_peak)

    if (plot == TRUE) {

      p <- histo %>%
        ggplot(aes(x = (n + c2_peak - 2*c1_peak) / (c2_peak - c1_peak))) +
        geom_line(aes(y = count), color = "red") +
        geom_line(aes(y = count_model), color = "blue") +
        # geom_line(aes(y = count_theo), color = "black") +
        # geom_vline(xintercept = 1:2, color = "black") +
        labs(x = "DNA content")

      if (!is.null(filename)) {
        ggsave(p, filename = filename)
      } else {
        p + labs(title = str_c( "B: ", B,
                                ", C: ", C,
                                ", D: ", tau - B - C,
                                ", tau: ", tau,
                                ", std_dev: ", std_dev,
                                ", std_dev_var: ", std_dev_var))
      }


    } else {
      return(tribble(~B, ~C, ~Std_dev, ~Std_dev_var, ~Deviation,
                     B,   C,  std_dev,  std_dev_var,  deviation))
    }
    
  }

### TEST 
# DNAhisto_currentGeneration(histo_test,
#                            B = 26, 
#                            C = 16, 
#                            tau = curr_tau, 
#                            std_dev = 17, 
#                            std_dev_var = 0.009, 
#                            chr_size = chr_size, 
#                            plot = TRUE)
## MODEL FITTING

### Loading histogram

histo_raw <- 
  histo_file %>% 
  read_tsv()


write(histo_file, stdout())

curr_reactor <- 
  histo_file %>% 
  str_match("\\/(.)_rep") %>% 
  .[2]

curr_tau <- 
  doubling_times %>% 
  dplyr::filter(reactor == curr_reactor) %>% 
  .$doubling_time

### Parameter grid search

percent_B <- seq(0, 1, 0.02)
percent_C <- seq(0, 1, 0.02)
std_dev_values <- seq(15, 19, 1)
std_dev_var_values <- seq(0.01, 0.05, 0.01)

### Fit

fits <- tribble(~B, ~C, ~Std_dev, ~Std_dev_var, ~Deviation)

fits_filename <- histo_file %>% str_replace("tsv", "fits")

i <- 0

for (b in percent_B) {
  
  for (c in percent_C) {
    
    if (b + c <= 1) {
      
      for (std_dev in std_dev_values) {
        
        for (std_dev_var in std_dev_var_values) {
          
          fits <- 
            DNAhisto_currentGeneration(histo_raw, 
                                       B = b * curr_tau, 
                                       C = c * curr_tau, 
                                       tau = curr_tau, 
                                       std_dev = std_dev, 
                                       std_dev_var = std_dev_var, 
                                       chr_size = chr_size) %>% 
            bind_rows(fits)
          
          i <- i + 1
          
          if (i %% 20 == 0) {
            fits %>% 
              arrange(Deviation) %>% 
              write_tsv(fits_filename, 
                        append = TRUE, 
                        col_names = !file.exists(fits_filename))
            
            fits <- tribble(~B, ~C, ~Std_dev, ~Std_dev_var, ~Deviation)
          }
          
        }
      }
    }
  }
}

### Best fit

# Reload all fits
# all_fits <- 
#   fits_filename %>% 
#   read_tsv(col_names = c("B", "C", "Std_dev", "Std_dev_var", "Deviation"))

all_fits <- 
  fits_file %>% 
  read_tsv(col_names = c("B", "C", "Std_dev", "Std_dev_var", "Deviation"),
          col_types = "ddidd")


# Plot best fit
best_fit <- 
  all_fits %>% 
  filter(B != 0, C != 0) %>% 
  top_n(n = -50, wt = Deviation)

best_fit_plotname <-
  histo_file %>% 
  str_replace("tsv", "best_fit.png")


DNAhisto_currentGeneration(histo_raw,
                           B = best_fit$B, 
                           C = best_fit$C, 
                           tau = curr_tau, 
                           std_dev = best_fit$Std_dev, 
                           std_dev_var = best_fit$Std_dev_var, 
                           chr_size = chr_size, 
                           plot = TRUE,
                           filename = best_fit_plotname)


# ## TEST

# DNAhisto_currentGeneration(histo_raw,
#                            B = 0.8 * curr_tau,
#                            C = 0.05 * curr_tau,
#                            tau = curr_tau,
#                            std_dev = 17,
#                            std_dev_var = 0.01,
#                            chr_size = chr_size,
#                            plot = TRUE)

c1 <- 18100
histo_raw %>% 
  ggplot(aes(pico, count)) + 
  geom_line() +
  geom_vline(xintercept = c1, color = "red") +
  geom_vline(xintercept = 2 * c1, color = "red")
  
