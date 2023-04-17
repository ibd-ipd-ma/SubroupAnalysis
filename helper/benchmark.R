# SOURCE: https://appsilon.com/fast-data-lookups-in-r-dplyr-vs-data-table/

time <- function(...) { 
  time_measurement <- system.time(eval(...)) 
  time_measurement[["user.self"]] 
} 

benchmark <- function(..., n = 10) { 
  times <- replicate(n, ...) 
  c(min = min(times), mean = mean(times), max = max(times)) 
}