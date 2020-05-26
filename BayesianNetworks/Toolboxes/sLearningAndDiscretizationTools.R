
# Local score functions ---------------------------------------------------

# Local score function, see Vogel(2013), p.33f (log deviation Monti(1998))

score_loc_bde_vogel = function(net_loc, dD_loc, dD_loc_i, n_iss, L_p){
  # Packages
  require('bnlearn')
  # Score
  score_loc = score(net_loc, dD_loc, 'bde', iss = n_iss) + # discrete part: Bayesian, dirichlet eqivalent uniform score
    L_p*drop(table(dD_loc_i)%*%log(1/table(dD_loc_i))) # continuous part
  return(score_loc)
}

score_loc_bic_vogel = function(net_loc, dD_loc, dD_loc_i, L_p){
  # Packages
  require('bnlearn')
  # Score
  score_loc = score(net_loc, dD_loc, 'bic') + # discrete part: Bayesian information criterion
    L_p*drop(table(dD_loc_i)%*%log(1/table(dD_loc_i))) # continuous part
  return(score_loc)
}

# Combining the two functions above in one function
score_loc_vogel = function(net_loc, dD_loc, dD_loc_i, score_name, n_iss, L_p){
  # Packages
  require('bnlearn')
  # Score
  if (score_name == 'bic' || score_name == 'bdj'){
    score_loc = score(net_loc, dD_loc, score_name) # discrete part: Bayesian, dirichlet eqivalent uniform score
  } else {
    score_loc = score(net_loc, dD_loc, score_name, iss = n_iss) # discrete part: Bayesian, dirichlet eqivalent uniform score
  }
  score_loc = score_loc + L_p*drop(table(dD_loc_i)%*%log(1/table(dD_loc_i))) # add continuous part 
  # Output
  return(score_loc)
}

# Local score function, see Monti(1998)

score_loc_monti = function(net_loc, dD_loc, dD_loc_i, D_loc_i, bp_list, score_name, n_iss, L_p){
  # Packages
  require('bnlearn')
  # Derived quantities
  mod_bp_list = bp_list;
  n_bp_list = length(mod_bp_list);
  mod_bp_list[1]=min(D_loc_i); mod_bp_list[n_bp_list]=max(D_loc_i); # avoid infinite interval length of boundery intervals
  int_lengths = mod_bp_list[-1] - mod_bp_list[-n_bp_list]; n_int = length(int_lengths);
  counts = table(dD_loc_i); n_count = length(counts);
  if (n_int != n_count){ # In case of two intervals with the same label, e.g. season variable [330;360]&[0;60] == "Winter"
    int_lengths[1] = int_lengths[1]+int_lengths[n_int] # combine the first and last interval
    int_lengths = int_lengths[-n_int]; # do.
  }
  # Score
  if (score_name == 'bic' || score_name == 'bdj'){
    score_loc = score(net_loc, dD_loc, score_name) # discrete part: Bayesian, dirichlet eqivalent uniform score
  } else {
    score_loc = score(net_loc, dD_loc, score_name, iss = n_iss) # discrete part: Bayesian, dirichlet eqivalent uniform score
  }
  score_loc = score_loc + L_p*drop(counts%*%log(1/int_lengths)) # add continuous part 
  # Output
  return(score_loc)
}


# Global score functions --------------------------------------------------

# Global score function, see Vogel(2013), p.33f (log deviation Monti(1998))

score_glo_vogel = function(net, d_data, score_name, n_iss, L_p){
  # Packages
  require('bnlearn')
  # Derived quantities
  nodes = names(d_data); # node names
  # Score
  if (score_name == 'bic' || score_name == 'bdj'){
    score_glo = score(net, d_data, score_name); # discrete part: Bayesian information criterion score
  } else {
    score_glo = score(net, d_data, score_name, iss = n_iss); # discrete part: Bayesian, dirichlet eqivalent uniform score
  }
  for (i in 1:length(nodes)){
    score_glo = score_glo + L_p*drop(table(d_data[nodes[i]])%*%log(1/table(d_data[nodes[i]]))); # add continuous part 
  }
  return(score_glo)
}

# Global score function, see Monti(1998)

score_glo_monti = function(net, d_data, data, bp_list, score_name, n_iss, L_p){
  # Packages
  require('bnlearn')
  # Derived quantities
  nodes = names(d_data); # node names
  # Score
  if (score_name == 'bic' || score_name == 'bdj'){
    score_glo = score(net, d_data, score_name); # discrete part: Bayesian information criterion score
  } else {
    score_glo = score(net, d_data, score_name, iss = n_iss); # discrete part: Bayesian, dirichlet eqivalent uniform score
  }
  for (i in 1:length(nodes)){
    mod_bp_list = bp_list[[nodes[i]]];
    n_bp_list = length(mod_bp_list);
    mod_bp_list[1]=min(data[nodes[i]]); mod_bp_list[n_bp_list]=max(data[nodes[i]]); # avoid infinite interval length of boundery intervals
    int_lengths = mod_bp_list[-1] - mod_bp_list[-n_bp_list]; n_int = length(int_lengths);
    counts = table(d_data[nodes[i]]); n_count = length(counts);
    if (n_int != n_count){ # In case two intervals with the same label, e.g. season variable [330;360]&[0;60] == "Winter"
      int_lengths[1] = int_lengths[1]+int_lengths[n_int] # combine the first and last interval
      int_lengths = int_lengths[-n_int]; # do.
    }
    # Updata score
    score_glo = score_glo + L_p*drop(counts%*%log(1/int_lengths)); # add continuous part 
  }
  return(score_glo)
}


# Discretization functions -------------------------------------------------

# My season discretizer ########################

# Discretizer 1
season_discretizer1 = function(vec, bp_list_bounds){
  # 'vec' is a vector in [0, 360] degrees indicating the season
  n_vec = length(vec)
  vec_fac = rep('0', n_vec);
  for (i in 1:n_vec){
    if ((vec[i] > 330) || (vec[i] <= 60)){
      vec_fac[i] = 'Winter'
    } else if ((vec[i] > 60) & (vec[i] <= 150)){
      vec_fac[i] = 'Spring'
    } else if ((vec[i] > 150) & (vec[i] <= 240)){
      vec_fac[i] = 'Summer'
    } else {
      vec_fac[i] = 'Fall'
    }
  }; vec_fac = factor(vec_fac) # store result as a factor
  bp_list = sort(c(bp_list_bounds, 60, 150, 240, 330))
  return(list(vec_fac = vec_fac, bp_list = bp_list))
}

# Discretizer 2
season_discretizer2 = function(vec, bp_list_bounds){
  # 'vec' is a vector in [0, 360] degrees indicating the season
  bp_list = sort(c(bp_list_bounds, 60, 150, 240, 330))
  vec_fac = cut(vec, breaks = bp_list)
  fac_levels = c('Winter', 'Spring', 'Summer', 'Fall', 'Winter')
  levels(vec_fac) = fac_levels;
  return(list(vec_fac = vec_fac, bp_list = bp_list))
}


# My direction discretizer ########################

direction_discretizer = function(vec, bp_list_bounds){
  # 'vec' is a vector in [0, 360] degrees indicating the direction
  bp_list = sort(c(bp_list_bounds, cumsum(c(22.5, rep(45,7)))))
  vec_fac = cut(vec, breaks = bp_list)
  fac_levels = c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N')
  levels(vec_fac) = fac_levels;
  return(list(vec_fac = vec_fac, bp_list = bp_list))
}


# My significant wave discretizer #################

Hm0_discretizer = function(vec, bp_list_bounds){
  # 'vec' is a vector of significant wave heights
  bp_list = sort(c(bp_list_bounds, cumsum(c(2, rep(0.5,14))))) # pr. 0.5 m (from 2 m to 9 m)
  vec_fac = cut(vec, breaks = bp_list)
  return(list(vec_fac = vec_fac, bp_list = bp_list))
}

Hm78_discretizer = function(vec, bp_list_bounds){
  # 'vec' is a vector of significant wave heights
  bp_list = sort(c(bp_list_bounds, cumsum(c(3, rep(1,13))))) # pr. 1 m (from 3 m to 16 m)
  vec_fac = cut(vec, breaks = bp_list)
  return(list(vec_fac = vec_fac, bp_list = bp_list))
}


# Discretization based on BP_list #################

BP_list_discretizer = function(data, d_data, nodes, bp_list){
  for (i in 1:length(nodes)){
    d_data[nodes[i]] = cut(data[nodes[i]][,1], breaks = bp_list[[nodes[i]]]) # discretize variable i
  }
  return(d_data) # save results
}


# My quantile discretizer ######################

quantile_discretizer = function(data, n_breaks, bp_list_bounds=NULL){
  # data: Data in a data.frame with named rows
  # n_breaks: number of intervals
  # bp_list_bounds: Upper and lower bounds for continuous variables, and bounds for descrite variables
  
  nodes = names(data)
  n_var = length(nodes)
  
  if (length(n_breaks) == 1){
    n_breaks = rep(n_breaks, ncol(data))
  }
  bp_list = list();
  d_data = data
  for (i in 1:n_var){
    name_i = nodes[i]
    data_i = data[name_i][,1]
    if (is(data_i, 'factor')){ # If a factor is included in the set of variables 
      stop("Variable ", nodes[i]," is a allready a factor. Remove this variable from the data set, and try again.")
    } else {
      quantiles = quantile((data_i), probs = seq(from = 0, to = n_breaks[i])/n_breaks[i])
      if (any(duplicated(quantiles))){
        stop("unable to discretize ", nodes[i], " in ", breaks[i],
             " intervals, some quantiles are not unique.")
      } # //*IF*//
      if (is.null(bp_list_bounds)==FALSE){
        quantiles[1] = unlist(bp_list_bounds[name_i], use.names = FALSE)[1] # replace lower bound
        quantiles[length(quantiles)] = unlist(bp_list_bounds[name_i], use.names = FALSE)[2] # replace upper bound  
      } # //*IF*//
      bp_list[[name_i]] = quantiles;
      d_data[[name_i]] = cut(data_i, breaks = quantiles)
    } # //*IF, ELSE*//
  } # //*FOR*//
  return(list(d_data=d_data, bp_list=bp_list))
} # //*FUNCTION*//

# Updata: handels partly observed data sets
quantile_discretizer_inclNA = function(data, n_breaks, bp_list_bounds=NULL){
  # data: Data in a data.frame with named rows
  # n_breaks: number of intervals
  # bp_list_bounds: Upper and lower bounds for continuous variables, and bounds for descrite variables
  
  nodes = names(data)
  n_var = length(nodes)
  
  if (length(n_breaks) == 1){
    n_breaks = rep(n_breaks, ncol(data))
  }
  bp_list = list();
  d_data = data
  for (i in 1:n_var){
    name_i = nodes[i]
    data_i = data[name_i][,1]
    if (is(data_i, 'factor')){ # If a factor is included in the set of variables 
      stop("Variable ", nodes[i]," is a allready a factor. Remove this variable from the data set, and try again.")
    } else {
      quantiles = quantile((data_i), probs = seq(from = 0, to = n_breaks[i])/n_breaks[i], na.rm = T)
      if (any(duplicated(quantiles))){
        stop("unable to discretize ", nodes[i], " in ", breaks[i],
             " intervals, some quantiles are not unique.")
      } # //*IF*//
      if (is.null(bp_list_bounds)==FALSE){
        quantiles[1] = unlist(bp_list_bounds[name_i], use.names = FALSE)[1] # replace lower bound
        quantiles[length(quantiles)] = unlist(bp_list_bounds[name_i], use.names = FALSE)[2] # replace upper bound  
      } # //*IF*//
      bp_list[[name_i]] = quantiles;
      d_data[[name_i]] = cut(data_i, breaks = quantiles)
    } # //*IF, ELSE*//
  } # //*FOR*//
  return(list(d_data=d_data, bp_list=bp_list))
} # //*FUNCTION*//


# Optimal discretization #####################

# One iteration over all variables ###
discretize_all_oneIt = function(BN, topo_order_con, DATA, dDATA, score_name_l, BP_list, BP_list_bounds, n_iss, L, 
                                quiet, n_iterations_max = NULL, score_con){
  # Packages
  require('bnlearn')
  
  n_var_con = length(topo_order_con) # no. of continuous variables to be discretized
  
  # Loop over all variables
  for (i in 1:n_var_con){ # Loop over all variables
    name_i = topo_order_con[i]; # name of variable i
    
    # Markov blanket of variable i
    MBincl_i = c(name_i , mb(BN, name_i)); # Variable 'i', and its Markov blanket (MB)
    local_net_i = subgraph(BN, MBincl_i) # Sub-graph where only variable 'i' and its MB is included
    
    dDATA_i = dDATA[MBincl_i] # Discrete data corr. to the variables in the local network
    data_i = DATA[name_i][,1] # Continuous data corr. to variable 'i'
    sort_data_i = sort(unique(data_i)) # sort, unique continuous data instances
    n_uniq_i = length(sort_data_i) # no. of unique continuous data instances
    BP_i = (sort_data_i[-n_uniq_i] + sort_data_i[-1])/2 # Initial set of potential boundery points (BP)
    if (length(BP_i)>512){
      BP_i = sort(sample(BP_i, 512, replace = FALSE, prob = NULL)); # Reduces the number of trial points
    }
    BP_list_i = unlist(BP_list_bounds[name_i], use.names = FALSE); # Initial list of BP 

    # Iteration for optimal discretization
    n_iterations = 0 # initialize number of iterations
    if (is.null(n_iterations_max)){ # maximum number of iterations (default) 
      ###################################
      if (round((n_uniq_i-1)/10) < 10){
        n_iterations_max_i = n_uniq_i-1 # min 10, but else ...
      } else {
        n_iterations_max_i = round((n_uniq_i-1)/10); # max 1/10 of potential BP
      }
      ###################################
    } else {
      ###################################
      if (length(n_iterations_max)==1){
        n_iterations_max_i = n_iterations_max
      } else {
        n_iterations_max_i = n_iterations_max[[name_i]];
      }
      ###################################
    } # //*IF, ELSE*//
    score_i = -Inf # Initialize reference score (iteration n_iter-1)
    score_i_new = -10^5 # Initialize new score
    
    while ((score_i - score_i_new < 0)&(n_iterations < n_iterations_max_i)) {
      n_iterations = n_iterations + 1; 
      
      if (quiet == FALSE){
        # print(c('Max_it',n_iterations_max_i)) 
        print(c('Var', 'Ite')); print(c(i,n_iterations)) # Number of iterations 
      } # //*IF*//
      
      # Updating
      if (n_iterations > 1){
        BP_i = BP_i_new
        BP_list_i = BP_list_i_new
        score_i = score_i_new 
        dDATA_i = dDATA_i_new 
      }
      
      # Loop over current set of potential BP
      score_ij_new = -Inf
      for (j in 1:length(BP_i)){ # Loop over boundery points (BP)
        
        BP_list_ij_test = sort(c(BP_list_i, BP_i[j])) # add new candidate point j to list of BP
        data_ij = cut(data_i, breaks = BP_list_ij_test) # discretize data according to new list
        dDATA_ij = dDATA_i; dDATA_ij[name_i] = data_ij; # Put new discrete variable in data set corr. to markov blanket
        
        # Score new configuration of BP's
        if (score_con == 'vogel') {
          score_ij = score_loc_vogel(net_loc = local_net_i,  dD_loc = dDATA_ij, dD_loc_i = data_ij, 
                                     score_name = score_name_l, n_iss = n_iss, L_p = L) # //[discretizer_tools.R]//
        } else if (score_con == 'monti') {
          score_ij = score_loc_monti(net_loc = local_net_i, dD_loc = dDATA_ij, dD_loc_i = data_ij, D_loc_i = data_i,
                                     bp_list = BP_list_ij_test, score_name = score_name_l, n_iss = n_iss, L_p = L) # //[discretizer_tools.R]//
        } # //*IF, ELSE IF*//
        
        # Update, iff score is improved compared to privious best candidate point
        if (score_ij > score_ij_new){ # Update, if the is an inprovenemt from privious trial
          index_ij_new = j # index of (current) best candidate point 
          BP_list_ij_new = BP_list_ij_test # new list of BP including this point
          score_ij_new = score_ij # new score
          dDATA_ij_new = dDATA_ij # new opdate date set
        } # //*IF*//
      } # //*FOR, j*//
      
      # Updating/book keeping
      BP_i_new = BP_i[-index_ij_new] # new list of candidate BP (exclude the one accepted in this iteration)  
      BP_list_i_new = BP_list_ij_new # new list of BP
      score_i_new = score_ij_new # new score
      dDATA_i_new = dDATA_ij_new # update data set
    }
    ###################################
    # Update discrete data
    if (n_iterations==n_iterations_max_i){
      dDATA[name_i] = dDATA_i_new[name_i]; # not cenverged  - choose last update! - ie. best possible score
      BP_list[[name_i]] = BP_list_i_new;
    } else{
      dDATA[name_i] = dDATA_i[name_i]; # converged - choose the minimum number of splits (one more reduces the score) - ie. best possible score
      BP_list[[name_i]] = BP_list_i;
    }
    ###################################
  } # //*FOR, i*//
  return(list(dDATA = dDATA, BP_list = BP_list))
} # //*FUNCTION*//

# Updata: handels partly observed data sets
discretize_all_oneIt_inclNA = function(BN, topo_order_con, DATA, DATAimp, dDATA, score_name_l, BP_list, BP_list_bounds, n_iss, L,
                                       quiet, n_iterations_max = NULL, score_con){
  # Packages
  require('bnlearn')
  
  n_var_con = length(topo_order_con) # no. of continuous variables to be discretized
  
  # Loop over all variables
  for (i in 1:n_var_con){ # Loop over all variables
    name_i = topo_order_con[i]; # name of variable i
    
    # Markov blanket of variable i
    MBincl_i = c(name_i , mb(BN, name_i)); # Variable 'i', and its Markov blanket (MB)
    local_net_i = subgraph(BN, MBincl_i) # Sub-graph where only variable 'i' and its MB is included
    
    dDATA_i = dDATA[MBincl_i] # Discrete data corr. to the variables in the local network
    data_i = DATA[name_i][,1] # Continuous data corr. to variable 'i' (incl. missings)
    data_imp_i = DATAimp[name_i][,1] # Continuous data corr. to variable 'i' (imputed)
    sort_data_i = sort(unique(na.omit(data_i))) # sort, unique continuous data instances (observed instances)
    n_uniq_i = length(sort_data_i) # no. of unique continuous data instances
    BP_i = (sort_data_i[-n_uniq_i] + sort_data_i[-1])/2 # Initial set of potential boundery points (BP)
    if (length(BP_i)>512){
      BP_i = sort(sample(BP_i, 512, replace = FALSE, prob = NULL)); # Reduces the number of trial points
    }
    BP_list_i = unlist(BP_list_bounds[name_i], use.names = FALSE); # Initial list of BP 
    
    # Iteration for optimal discretization
    n_iterations = 0 # initialize number of iterations
    if (is.null(n_iterations_max)){ # maximum number of iterations (default) 
      ###################################
      if (round((n_uniq_i-1)/10) < 10){
        n_iterations_max_i = n_uniq_i-1 # 1/10 of potential BP is less than 10
      } else {
        n_iterations_max_i = round((n_uniq_i-1)/10); # max 1/10 of potential BP
      } 
      ###################################
    } else {
      ###################################
      if (length(n_iterations_max)==1){
        n_iterations_max_i = n_iterations_max
      } else {
        n_iterations_max_i = n_iterations_max[[name_i]];
      }
      ###################################
    } # //*IF, ELSE*//
    score_i = -Inf # Initialize reference score (iteration n_iter-1)
    score_i_new = -10^5 # Initialize new score
    
    while ((score_i - score_i_new < 0)&(n_iterations < n_iterations_max_i)) {
      n_iterations = n_iterations + 1; 
      
      if (quiet == FALSE){
        # print(c('Max_it',n_iterations_max_i)) 
        print(c('Var', 'Ite')); print(c(i,n_iterations)) # Number of iterations 
      } # //*IF*//
      
      # Updating
      if (n_iterations > 1){
        BP_i = BP_i_new
        BP_list_i = BP_list_i_new
        score_i = score_i_new 
        dDATA_i = dDATA_i_new 
      }
      
      # Loop over current set of potential BP
      score_ij_new = -Inf
      for (j in 1:length(BP_i)){ # Loop over boundery points (BP)
        
        BP_list_ij_test = sort(c(BP_list_i, BP_i[j])) # add new candidate point j to list of BP
        data_ij = cut(data_imp_i, breaks = BP_list_ij_test) # discretize data according to new list
        dDATA_ij = dDATA_i; dDATA_ij[name_i] = data_ij; # Put new discrete variable in data set corr. to markov blanket
        
        # Score new configuration of BP's
        if (score_con == 'vogel') {
          score_ij = score_loc_vogel(net_loc = local_net_i,  dD_loc = dDATA_ij, dD_loc_i = data_ij, 
                                     score_name = score_name_l, n_iss = n_iss, L_p = L) # //[discretizer_tools.R]//
        } else if (score_con == 'monti') {
          score_ij = score_loc_monti(net_loc = local_net_i, dD_loc = dDATA_ij, dD_loc_i = data_ij, D_loc_i = data_imp_i,
                                     bp_list = BP_list_ij_test, score_name = score_name_l, n_iss = n_iss, L_p = L) # //[discretizer_tools.R]//
        } # //*IF, ELSE IF*//
        
        # Update, iff score is improved compared to privious best candidate point
        if (score_ij > score_ij_new){ # Update, if the is an inprovenemt from privious trial
          index_ij_new = j # index of (current) best candidate point 
          BP_list_ij_new = BP_list_ij_test # new list of BP including this point
          score_ij_new = score_ij # new score
          dDATA_ij_new = dDATA_ij # new opdate date set
        } # //*IF*//
      } # //*FOR, j*//
      
      # Updating/book keeping
      BP_i_new = BP_i[-index_ij_new] # new list of candidate BP (exclude the one accepted in this iteration)  
      BP_list_i_new = BP_list_ij_new # new list of BP
      score_i_new = score_ij_new # new score
      dDATA_i_new = dDATA_ij_new # update data set
    }
    ###################################
    # Update discrete data
    if (n_iterations==n_iterations_max_i){
      dDATA[name_i] = dDATA_i_new[name_i]; # not cenverged  - choose last update! - ie. best possible score
      BP_list[[name_i]] = BP_list_i_new;
    } else{
      dDATA[name_i] = dDATA_i[name_i]; # converged - choose the minimum number of splits (one more reduces the score) - ie. best possible score
      BP_list[[name_i]] = BP_list_i;
    }
    ###################################
  } # //*FOR, i*//
  return(list(dDATA = dDATA, BP_list = BP_list))
} # //*FUNCTION*//


# Imputation --------------------------------------------------------------

# Imputation of continuous values based on discrete levels ############
impute_continuous = function(data, ddata, bp_list){
  data_imp = data;
  nodes = names(data);
  for (i in 1:length(nodes)){ # loop over variables
    name_i = nodes[i]
    data_i = data[name_i]
    if (anyNA(data_i)==TRUE){ # test, if variable contain NAs 
      na_index_i = which(is.na(data_i)) # index of missings
      bp_list_i = unlist(bp_list[name_i], use.names = FALSE); # Define boundery points
      bp_list_i[1] = min(data_i, na.rm = T); bp_list_i[length(bp_list_i)] = max(data_i, na.rm = T); # Makes sure that continuous imputations are not unbounded
      ddata_i = ddata[name_i]; # imputed, discrete data set
      
      int_levels_i = as.integer(ddata_i[na_index_i,1]) # find discrete values of missings
      uniq_levels_i = sort(unique(int_levels_i)) # define unique levels
      
      for (j in 1:length(uniq_levels_i)){ # loop over unique levels
        index_level_ij = which(int_levels_i==uniq_levels_i[j]) # difine missings belonging to class, j
        n_levels_ij = length(index_level_ij) # no. of missings belonging to class, j
        data_ij = runif(n_levels_ij, min = bp_list_i[uniq_levels_i[j]], max = bp_list_i[uniq_levels_i[j]+1]) # impute missings labeled, j (based on generative model)
        data_i[na_index_i[index_level_ij],1] = data_ij # replace missing values in varaible, i
      }
    }
    data_imp[name_i] = data_i # define an imputed data set
  }
  return(data_imp)
}

# Updata: Includes bootstrap resampling for continious values
impute_continuous_v2 = function(data, ddata, bp_list, score_con){
  data_imp = data;
  nodes = names(data);
  for (i in 1:length(nodes)){ # loop over variables
    name_i = nodes[i]
    data_i = data[name_i]
    if (anyNA(data_i)==TRUE){ # test, if variable contain NAs 
      na_index_i = which(is.na(data_i)) # index of missings
      bp_list_i = unlist(bp_list[name_i], use.names = FALSE); # Define boundery points
      bp_list_i[1] = min(data_i, na.rm = T); bp_list_i[length(bp_list_i)] = max(data_i, na.rm = T); # Makes sure that continuous imputations are not unbounded
      ddata_i = ddata[name_i]; # imputed, discrete data set
      
      int_levels_i = as.integer(ddata_i[na_index_i,1]) # find discrete values of missings
      uniq_levels_i = sort(unique(int_levels_i)) # define unique levels
      
      data_i_obs = data_i[-na_index_i]
      for (j in 1:length(uniq_levels_i)){ # loop over unique levels
        index_level_ij = which(int_levels_i==uniq_levels_i[j]) # difine missings belonging to class, j
        n_levels_ij = length(index_level_ij) # no. of missings belonging to class, j
        
        data_i_obs = data_i[-na_index_i, ]
        #################################
        if (score_con=='monti'){
          data_ij = runif(n_levels_ij, min = bp_list_i[uniq_levels_i[j]], max = bp_list_i[uniq_levels_i[j]+1]) # impute missings labeled, j (based on generative model)  
        } else {
          samp_ij = data_i_obs[ (data_i_obs >= bp_list_i[uniq_levels_i[j]]) & (data_i_obs <= bp_list_i[uniq_levels_i[j]+1]) ]
          data_ij = sample(x=samp_ij, size = n_levels_ij, replace = T) # impute missings labeled, j (based on generative model)  
        }
        #################################
        
        data_i[na_index_i[index_level_ij],1] = data_ij # replace missing values in varaible, i
      }
    }
    data_imp[name_i] = data_i # define an imputed data set
  }
  return(data_imp)
}


# Structure Learning ------------------------------------------------------

# Learn structure and MAP parameters, as well as most likely imputed data set
LearnStructure_inclNA = function(dDATA, search_alg, score_name_g, n_iss, WL, BL, imp_method, net_start){
  # Packages
  require('bnlearn')
  
  # Structure learning
  if (search_alg == 'hill'){
    if (score_name_g == 'bic' || score_name_g == 'bdj'){
      RESstruc = structural.em(dDATA, maximize = 'hc', maximize.args = list(score = score_name_g, restart = 5, perturb = 5, whitelist = WL, blacklist = BL), 
                               fit = "bayes", fit.args = list(iss = n_iss), impute = imp_method, start = net_start, return.all = TRUE) 
    } else { # 'bde'
      RESstruc = structural.em(dDATA, maximize = 'hc', maximize.args = list(score = score_name_g, restart = 5, perturb = 5, iss = n_iss, whitelist = WL, blacklist = BL), 
                               fit = "bayes", fit.args = list(iss = n_iss), impute = imp_method, start = net_start, return.all = TRUE)
    } # //*IF, ELSE IF; 'score_name_g'*//
  } else if (search_alg == 'tabu'){
    if (score_name_g == 'bic' || score_name_g == 'bdj'){
      RESstruc = structural.em(dDATA, maximize = 'tabu', maximize.args = list(score = score_name_g, whitelist = WL, blacklist = BL), 
                               fit = "bayes", fit.args = list(iss = n_iss), impute = imp_method, start = net_start, return.all = TRUE)
    } else { # 'bde'
      RESstruc = structural.em(dDATA, maximize = 'tabu', maximize.args = list(score = score_name_g, iss = n_iss, whitelist = WL, blacklist = BL), 
                               fit = "bayes", fit.args = list(iss = n_iss), impute = imp_method, start = net_start, return.all = TRUE)
    } # //*IF, ELSE IF; 'score_name_g'*//
  } # //*IF, ELSE IF; 'search_alg'*//
  return(RESstruc)
}


# Combined discretization and learning  -----------------------------------

discretizeANDlearn_all = function(score_name_g, score_name_l, L, quiet, score_con, search_alg,
                                  DATA, dDATA, BP_list, BP_list_bounds, n_iterations_loc_max = NULL, 
                                  var_con, net_start, WL = NULL, BL = NULL, scale_iss){
  # Packages
  require('bnlearn')
  
  nodes = names(DATA) # variable names
  
  if (quiet == FALSE){
    if (is.null(net_start)){
      plot(empty.graph(nodes)) # plot network used to start the initial search
    } else {
      plot(net_start) # plot network used to start the initial search
    }
  }
  
  # Initial network
  if (search_alg == 'hill'){
    if (score_name_g == 'bic' || score_name_g == 'bdj'){
      n_iss = NULL; rm(scale_iss)
      iBN = hc(dDATA, score = score_name_g, start = net_start, whitelist = WL, blacklist = BL, restart = 5, perturb = 5) # initial BN found by hill-climbing, and bayesian information criterion
    } else { # 'bde'
      n_iss = 1*scale_iss; # imaginary sample size in dirichlet distribution
      iBN = hc(dDATA, score = score_name_g, iss = n_iss, start = net_start, whitelist = WL, blacklist = BL, restart = 5, perturb = 5) # initial BN found by hill-climbing, and bayesian, dirichlet equivalent uniform score 
    } # //*IF, ELSE IF; 'score_name_g'*//
  } else if (search_alg == 'tabu'){
    if (score_name_g == 'bic' || score_name_g == 'bdj'){
      n_iss = NULL; rm(scale_iss)
      iBN = tabu(dDATA, score = score_name_g, start = net_start, whitelist = WL, blacklist = BL) # initial BN found by tabu-search, and bayesian information criterion
    } else { # 'bde'
      n_iss = 1*scale_iss; # imaginary sample size in dirichlet distribution
      iBN = tabu(dDATA, score = score_name_g, iss = n_iss, start = net_start, whitelist = WL, blacklist = BL) # initial BN found by tabu_search, and bayesian, dirichlet equivalent uniform score 
    } # //*IF, ELSE IF; 'score_name_g'*//
  } # //*IF, ELSE IF; 'search_alg'*//
  if (quiet == FALSE){
    plot(iBN) # plot initial network
  }
  
  # Discretizer #######################
  
  # Initialization
  BN = iBN # initial network
  topo_order = node.ordering(BN) # initial topological ordering
  topo_order_con = intersect(topo_order, var_con) # Continuous variables in topological order
  
  if (score_con == 'vogel'){
    score_glo = score_glo_vogel(net = BN, d_data = dDATA, score_name = score_name_g, # Initial score of BN and discretization
                                n_iss = n_iss, L_p = L) # //*DISCRETIZER_TOOLS.R*//
  } else if (score_con == 'monti'){
    score_glo = score_glo_monti(net = BN, d_data = dDATA, data = DATA, bp_list = BP_list, # Initial score of BN and discretization
                                score_name = score_name_g, n_iss = n_iss, L_p = L) # //*DISCRETIZER_TOOLS.R*//
  }
  
  n_iterations_glo = 0; # initialize no. of iterations
  n_iterations_glo_max = length(nodes); # max no. of iterations
  score_glo_new = 0; # initialize new global score
  while((n_iterations_glo <= n_iterations_glo_max)&(abs(score_glo_new/score_glo - score_glo/score_glo_new) > 0.0010)){ # Convergence
    n_iterations_glo = n_iterations_glo + 1; # update iteration no.
    
    if (quiet == FALSE){
      print(c('Global iteration', n_iterations_glo)) # print progress
    }
    
    # Update (after first run)
    if (n_iterations_glo > 1){
      score_glo = score_glo_new; # update baseline score
      dDATA = dDATA_new; # update discrete, baseline data set
      BP_list = BP_list_new; # update list of boundery points
      BN = BN_new; # update baseline BN
      topo_order = node.ordering(BN) # update topological ordering
      topo_order_con = intersect(topo_order, var_con) # Continuous variables in topological order
    }
    
    # Re-discretize continuous variables
    res_oneIt = discretize_all_oneIt(BN = BN, topo_order_con = topo_order_con, DATA = DATA, 
                                     dDATA = dDATA, score_name_l = score_name_l, BP_list = BP_list, 
                                     BP_list_bounds = BP_list_bounds, n_iss = n_iss, L = L, quiet = quiet, 
                                     n_iterations_max = n_iterations_loc_max, score_con = score_con) # //*DISCRETIZER_TOOLS.R*//
    # Update variables
    dDATA_new = res_oneIt$dDATA;
    BP_list_new = res_oneIt$BP_list;
    
    # New network
    if (search_alg == 'hill'){
      if (score_name_g == 'bic' || score_name_g == 'bdj'){
        BN_new = hc(dDATA_new, score = score_name_g, start = BN, whitelist = WL, blacklist = BL, restart = 5, perturb = 5)
      } else { # 'bde'
        BN_new = hc(dDATA_new, score = score_name_g, iss = n_iss, start = BN, whitelist = WL, blacklist = BL, restart = 5, perturb = 5) 
      } # // IF, ELSE IF; 'score_name_g'*//
    } else if (search_alg == 'tabu'){
      if (score_name_g == 'bic' || score_name_g == 'bdj'){
        BN_new = tabu(dDATA_new, score = score_name_g, start = BN, whitelist = WL, blacklist = BL)
      } else { # 'bde'
        BN_new = tabu(dDATA_new, score = score_name_g, iss = n_iss, start = BN, whitelist = WL, blacklist = BL) 
      } # // IF, ELSE IF; 'score_name_g'*//
    } # // IF, ELSE IF; 'search_alg'*//
    
    # New network score
    if (score_con == 'vogel'){
      score_glo_new = score_glo_vogel(net = BN_new, d_data = dDATA_new, score_name = score_name_g, 
                                      n_iss = n_iss, L_p = L) # //*DISCRETIZER_TOOLS.R*//
    } else if (score_con == 'monti'){
      score_glo_new = score_glo_monti(net = BN_new, d_data = dDATA_new, data = DATA, bp_list = BP_list_new, # Initial score of BN and discretization
                                      score_name = score_name_g, n_iss = n_iss, L_p = L) # //*DISCRETIZER_TOOLS.R*//
    }
    # Progress
    if (quiet == FALSE){
      plot(BN_new) # plot proposal network
      print(c('old score', score_glo)) # print score
      print(c('new score', score_glo_new)) # print new score
    }
  }
  return(list(BN = BN_new, dDATA = dDATA_new, BP_list = BP_list_new))
}


# Updata: handels partly observed data sets
discretizeANDlearn_all_inclNA = function(score_name_g, score_name_l, L, quiet, score_con, search_alg, imp_method,
                                         DATA, dDATA, BP_list, BP_list_bounds, n_iterations_loc_max = NULL, var_con, 
                                         net_start, WL = NULL, BL = NULL, n_iss){
  # Packages
  require('bnlearn')
  
  nodes = names(DATA) # variable names
  id_missing = arrayInd(which(is.na(DATA)), dim(DATA)) # indeices of missing data
  
  # Initial network
  RESstruc = LearnStructure_inclNA(dDATA=dDATA, search_alg=search_alg, score_name_g=score_name_g, n_iss=n_iss, WL=WL, BL=BL, 
                                   imp_method = imp_method, net_start=net_start) # //*DISCRETIZER_TOOLS.R*//
  if (quiet == FALSE){
    plot(RESstruc$dag) # plot initial network
  }
  
  # Discretizer #######################
  
  # Initialization
  dDATAimp = RESstruc$imputed # completed (discrete) data set
  ###################################
  # DATAimp = impute_continuous(data = DATA, ddata = dDATAimp, bp_list = BP_list) # completed (continuous) data set
  DATAimp = impute_continuous_v2(data = DATA, ddata = dDATAimp, bp_list = BP_list, score_con = score_con)
  ###################################
  BN = RESstruc$dag # initial network
  topo_order = node.ordering(BN) # initial topological ordering
  topo_order_con = intersect(topo_order, var_con) # Continuous variables in topological order
  
  if (score_con == 'vogel'){
    score_glo = score_glo_vogel(net = BN, d_data = dDATAimp, score_name = score_name_g, # Initial score of BN and discretization
                                n_iss = n_iss, L_p = L) # //*DISCRETIZER_TOOLS.R*//
  } else if (score_con == 'monti'){
    score_glo = score_glo_monti(net = BN, d_data = dDATAimp, data = DATAimp, bp_list = BP_list, # Initial score of BN and discretization
                                score_name = score_name_g, n_iss = n_iss, L_p = L) # //*DISCRETIZER_TOOLS.R*//
  }
  
  # Iteration procedure
  n_iterations_glo = 0; # initialize no. of iterations
  n_iterations_glo_max = length(nodes); # max no. of iterations
  score_glo_new = 0; # initialize new global score
  while((n_iterations_glo <= n_iterations_glo_max)&(abs(score_glo_new/score_glo - score_glo/score_glo_new) > 0.0010)){ # Convergence
    n_iterations_glo = n_iterations_glo + 1; # update iteration no.
    
    if (quiet == FALSE){
      print(c('Global iteration', n_iterations_glo)) # print progress
    }
    
    # Update (after first run)
    if (n_iterations_glo > 1){
      score_glo = score_glo_new; # update baseline score
      dDATA = dDATA_new; # update discrete, baseline data set
      DATAimp = DATAimp_new;
      dDATAimp = dDATAimp_new;
      BP_list = BP_list_new; # update list of boundery points
      BN = BN_new; # update baseline BN
      topo_order = node.ordering(BN) # update topological ordering
      topo_order_con = intersect(topo_order, var_con) # Continuous variables in topological order
    }
    
    # Re-discretize continuous variables
    res_oneIt = discretize_all_oneIt_inclNA(BN = BN, topo_order_con = topo_order_con, DATA = DATA, DATAimp = DATAimp,
                                            dDATA = dDATAimp, score_name_l = score_name_l, BP_list = BP_list,
                                            BP_list_bounds = BP_list_bounds, n_iss = n_iss, L = L, quiet = quiet, 
                                            n_iterations_max = n_iterations_loc_max, score_con = score_con) # //*DISCRETIZER_TOOLS.R*//
    # Update variables
    dDATA_new = res_oneIt$dDATA;
    dDATA_new[id_missing] = NA;
    BP_list_new = res_oneIt$BP_list;
    
    # Updata network
    BNfit_start = bn.fit(BN, res_oneIt$dDATA, method = 'bayes', iss = n_iss)
    RESstruc_new = LearnStructure_inclNA(dDATA=dDATA_new, search_alg=search_alg, score_name_g=score_name_g, n_iss=n_iss, WL=WL, BL=BL, 
                                         imp_method = imp_method, net_start=BNfit_start) # //*DISCRETIZER_TOOLS.R*//
    
    dDATAimp_new = RESstruc_new$imputed # completed (discrete) data set
    ############################
    # DATAimp_new = impute_continuous(data = DATA, ddata = dDATAimp_new, bp_list = BP_list_new) # completed (continuous) data set
    DATAimp_new = impute_continuous_v2(data = DATA, ddata = dDATAimp_new, bp_list = BP_list_new, score_con = score_con) # completed (continuous) data set
    ############################
    BN_new = RESstruc_new$dag # initial network
    
    # New network score
    if (score_con == 'vogel'){
      score_glo_new = score_glo_vogel(net = BN_new, d_data = dDATAimp_new, score_name = score_name_g, 
                                      n_iss = n_iss, L_p = L) # //*DISCRETIZER_TOOLS.R*//
    } else if (score_con == 'monti'){
      score_glo_new = score_glo_monti(net = BN_new, d_data = dDATAimp_new, data = DATAimp_new, bp_list = BP_list_new, # Initial score of BN and discretization
                                      score_name = score_name_g, n_iss = n_iss, L_p = L) # //*DISCRETIZER_TOOLS.R*//
    }
    # Progress
    if (quiet == FALSE){
      plot(BN_new) # plot proposal network
      print(c('old score', score_glo)) # print score
      print(c('new score', score_glo_new)) # print new score
    }
  }
  return(list(BN = BN_new, dDATA = dDATA_new, BP_list = BP_list_new, 
              dDATAimp = dDATAimp_new, DATAimp = DATAimp_new, FITTED = RESstruc_new$fitted))
}