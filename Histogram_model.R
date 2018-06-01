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
             "Usage: Rscript --vanilla run_histogram_model.R <full input file path>"),
       call. = FALSE)
}

histo_file = args[1]

## CONSTANTS

chr_size <- 1.3

doubling_times <- #doubling times in days
  tibble(reactor = c("A", "Z", "M", "R"),
         doubling_time = c(12.1, 38.8, 48.5, 23.1))


## FUNCTIONS 

### Get main peak channel

getC1 <- function(data_binned) {
  # Given a df of binned DNA data (channel ~ count), returns the C1 channel
  c1_bin <- 
    data_binned %>% 
    top_n(5, count) %>% 
    .$n %>% 
    mean() %>% 
    round()
  
  return(c1_bin)
}


### DNA theoretical distribution

#### Theoretical count in a given channel

theoretical_count <- function(bin, start_bin, B, C, tau, total_n, time_per_channel) {
  # Given a c1 channel, B, C, tau and total counts, returns the theoretical counts
  # in a specific channel.
  
  if (bin < start_bin | bin > 2 * start_bin) {
    count_theo <- 0
  } else if (bin == start_bin) {
    count_theo <- total_n * pop_percent(0, B, tau)
  } else if (bin == 2 * start_bin) {
    count_theo <- total_n * pop_percent(B + C, tau, tau)
  } else {
    count_theo <- total_n * pop_percent(B + time_per_channel * (bin - start_bin),
                                        B + time_per_channel * (bin + 1 - start_bin),
                                        tau) # might be missing a bin during intergration
  }
  
  return(count_theo)
}

#### Percentage of cells in a given time slice

pop_percent <- function(t1, t2, t_total) {
  2 * (2 ^ (-t1 / t_total) - 2 ^ (-t2 / t_total))
}


### DNA histogram model for replication in current generation (slow growth)

DNAhisto_currentGeneration <- 
  function(histo, B, C, tau, std_dev, std_dev_var, chr_size, plot = FALSE, filename = "") {
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
    
    if ( B + C > tau) {
      print(str_c("B: ", B))
      print(str_c("C: ", C))
      print("DNA replication didn't start in current generation.")
      return(tribble(~B, ~C, ~Std_dev, ~Std_dev_var, ~Deviation,
                      B,  C,  std_dev,  std_dev_var,  NA))
    }
    
    c1_peak <- getC1(histo) 
    total_count <- histo$count %>% sum()
    dna_per_bin <- chr_size / c1_peak
    start_bin <- chr_size / dna_per_bin
    time_per_channel <- C / start_bin
    
    
    histo <-
      histo %>%
      mutate(count_theo = map_dbl(n, 
                                  ~theoretical_count(., start_bin, B, C, tau, total_count, time_per_channel))
      )
    
    # std_dev with incremental change
    histo <-
      histo %>% 
      mutate(std_dev = std_dev + std_dev_var * row_number(),
             count_model = map2(n, std_dev, 
                                function(x, y) {histo$count_theo * dnorm(histo$n, x, y) %>% 
                                    as_tibble}) %>% 
               bind_cols() %>% 
               summarise_all(sum) %>% 
               gather(key = x, value = y) %>% 
               .$y)
    
    deviation <- 
      histo %>% 
      mutate(dif = (sqrt(count) - sqrt(count_model))^2) %>% 
      # filter(n > 0.5 * c1_peak) %>% # ignore the debris / cells without full DNA
      transmute(dif = dif) %>% 
      sum() %>% 
      sqrt(.) / 255
    
    if (plot == TRUE) {
      
      p <- histo %>% 
        ggplot(aes(x = n / c1_peak)) +
        geom_line(aes(y = count), color = "red") +
        geom_line(aes(y = count_model), color = "blue") +
        annotate("text", x = 2.3, y = max(0.9 * histo$count), 
                 label = stringr::str_c("Cell cycle (in days)\nB: ", B,
                                        "\nC: ", formatC(C, digits = 2, format = "f"),
                                        "\nD: ", formatC((tau - C - B), digits = 2, format = "f"),
                                        "\ndeviation: ", formatC(deviation, digits = 2, format = "f"))) +
        labs(x = "DNA content")
      
      ggsave(p, filename = filename)
      
    } else {
      return(tribble(~B, ~C, ~Std_dev, ~Std_dev_var, ~Deviation,
                     B,   C,  std_dev,  std_dev_var,  deviation))
    }
    
  }

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

tau <- 
  doubling_times %>% 
  dplyr::filter(reactor == curr_reactor) %>% 
  .$doubling_time

### Parameter grid search

percent_B <- seq(0, 1, 0.2)
percent_C <- seq(0, 1, 0.2)
std_dev_values <- seq(5, 20, 5)
std_dev_var_values <- seq(0.0, 0.15, 0.05)

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
                                       B = b * tau, 
                                       C = c * tau, 
                                       tau = tau, 
                                       std_dev = std_dev, 
                                       std_dev_var = std_dev_var, 
                                       chr_size = chr_size) %>% 
            bind_rows(fits)
          
          i <- i + 1
          
          if (i %% 10 == 0) {
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
all_fits <- 
  fits_filename %>% 
  read_tsv(col_names = c("B", "C", "Std_dev", "Std_dev_var", "Deviation"))


# Plot best fit
best_fit <- 
  all_fits %>% 
  top_n(n = -1, wt = Deviation)

best_fit_plotname <-
  histo_file %>% 
  str_replace("tsv", "best_fit.png")


DNAhisto_currentGeneration(histo_raw,
                           B = best_fit$B, 
                           C = best_fit$C, 
                           tau = tau, 
                           std_dev = best_fit$Std_dev, 
                           std_dev_var = best_fit$Std_dev_var, 
                           chr_size = chr_size, 
                           plot = TRUE,
                           filename = best_fit_plotname)
