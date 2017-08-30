
# A HCR is precautionary if the probability that B < B_lim is less than 5%.
# Other performance criteria: long-term average yield, variance of yield, 
# time to recovery to some biomass level



is_precautionary(MSE) {
  if(!inherits(MSE, "MSE")) stop("No object of class MSE was detected")
  
  return(1)
}