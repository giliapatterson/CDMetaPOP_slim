grow_temperature <- function(current_size, age_class, patch_temp, grow_days, popvars){
  X_max = popvars$growth_temp_max
  X_CV = popvars$growth_temp_CV
  R0 = popvars$growth_R0
  Loo = popvars$growth_Loo
  t0 = popvars$growth_temp_t0
  
  int_R = -R0*dnorm(patch_temp, X_max, X_CV*X_max)/dnorm(X_max,  X_max, X_CV*X_max);
  L_inc = Loo * (1-exp(int_R*(age_class + 1-t0)))* dnorm(patch_temp, X_max, X_CV*X_max)/dnorm(X_max,  X_max, X_CV*X_max);
  L_inc_age = L_inc*exp((age_class+1)*(-R0));
  newsize = current_size + (L_inc_age * (grow_days/365));
  return(newsize);
}

# Calculate probability of maturity
prob_mature <- function(size, sex_string, age, popvars){
  params <- popvars |>
    mutate(sex = "F~M") |>
    separate_longer_delim(everything(), "~") |>
    filter(sex == sex_string) |>
    mutate(across(c(mature_age, mature_eqn_int, mature_eqn_slope), as.numeric))
  if(age >= params$mature_age){
    p_mature = 1.0
  }
  else{
    p_mature = exp(params$mature_eqn_int + params$mature_eqn_slope*size) / (1 + exp(params$mature_eqn_int + params$mature_eqn_slope*size))
  }
  return(p_mature);
}

prob_mature_vector <- function(sizes, sexes, ages, popvars){
  prob_mature_res <- rep(0, length(sizes))
  for(i in 1:length(sizes)){
    if(length(sexes) == 1){
      sex = sexes
    }   
    if(length(sexes) > 1){
      sex = sexes[i]
    }   
    prob_mature_res[i] = prob_mature(sizes[i], sex, ages[i], popvars)
  }
  return(prob_mature_res)
}

# Number of eggs produced by an individual
mean_num_eggs <- function(sizes, popvars){
  popvars <- popvars |>
    mutate(across(c(Egg_Mean_par1, Egg_Mean_par2) , as.numeric))
  E_mu = rep(0, length(sizes))
  for(i in 1:length(sizes)){
    E_mu[i] = popvars$Egg_Mean_par1*exp(popvars$Egg_Mean_par2)*sizes[i]
    #noffspring = rpois(1, E_mu);
  }
  return(E_mu);
}