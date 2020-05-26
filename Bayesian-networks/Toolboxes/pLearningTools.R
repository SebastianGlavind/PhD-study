# Bayesian Dirichlet equivalent -------------------------------------------

# BDeu prior (not used below; instead it is embedded)
BDeu_prior = function(BN, dDATA, n_iss){
  require(bnlearn)
  
  nodes = names(dDATA)
  PCprior=list()
  for (i in 1:length(nodes)){
    nodes_i = nodes[i]
    parents_i = bnlearn::parents(BN, nodes_i); n_pa_i = length(parents_i)
    levels_nodes_i = levels(dDATA[nodes_i][1,]);
    dim_iCT = length(levels_nodes_i)
    dim_iNames = list(); dim_iNames[[nodes_i]] = levels_nodes_i
    if (n_pa_i!=0){
      for (j in 1:n_pa_i){
        levels_pa_ij = levels(dDATA[parents_i[j]][1,])
        dim_iCT = c(dim_iCT, length(levels_pa_ij))
        dim_iNames[[parents_i[j]]] = levels_pa_ij
      }
    }
    n_iEnteries = prod(dim_iCT); # enteries in iCT
    iPseC_scale = n_iss/n_iEnteries; # Bayesian Dirichlet equivalent (BDe); prior scale
    iPseC = as.table(array(iPseC_scale, dim=dim_iCT, dimnames = dim_iNames)) # BDe prior counts (in table form)
    PCprior[[nodes_i]] = iPseC
  }
  return(PCprior)
}

# Posterior counts - input is bnfit object
posterior_counts = function(bn_map, ddata, n_iss){
  PC=list();
  for (i in 1:length(bn_map)){
    iVar = bn_map[[i]]
    iNode = iVar$node;
    iPar = iVar$parents; n_iPar=length(iPar);
    
    # Contingency table (data)
    iddata = ddata[iNode];
    if (n_iPar != 0){
      iddata = cbind(iddata,  ddata[iPar])
      iCT=table(iddata) # contingency table
    } else {
      iCT=table(iddata, dnn = iNode) # contingency table
      # iCT=table(iddata, deparse.level = 0) # contingency table
    }
    
    # Pseudo counts (imaginary sample)
    dim_iCT = dim(iCT) # dimention of iCT
    n_iEnteries = prod(dim_iCT); # enteries in iCT
    iPseC_scale = n_iss/n_iEnteries; # Bayesian Dirichlet equivalent (BDe); prior scale
    iPseC = array(iPseC_scale,dim_iCT) # BDe prior counts (in table form)
    
    # BDe; posterior counts
    iPC = iCT + iPseC 
    
    # store results
    PC[[iNode]] = iPC;
  }
  return(PC)
}

# Posterior counts - input is bn object
posterior_counts_bn = function(bn, ddata, n_iss){
  require(bnlearn)
  
  PC=list();
  nodes = names(ddata)
  for (i in 1:length(nodes)){
    iNode = nodes[i]
    iPar = bnlearn::parents(bn,  iNode); n_iPar=length(iPar);
    
    # Contingency table (data)
    iddata = ddata[iNode];
    if (n_iPar != 0){
      iddata = cbind(iddata,  ddata[iPar])
      iCT=table(iddata) # contingency table
    } else {
      iCT=table(iddata, dnn = iNode) # contingency table
      # iCT=table(iddata, deparse.level = 0) # contingency table
    }
    
    # Pseudo counts (imaginary sample)
    dim_iCT = dim(iCT) # dimention of iCT
    n_iEnteries = prod(dim_iCT); # enteries in iCT
    iPseC_scale = n_iss/n_iEnteries; # Bayesian Dirichlet equivalent (BDe); prior scale
    iPseC = array(iPseC_scale,dim_iCT) # BDe prior counts (in table form)
    
    # BDe; posterior counts
    iPC = iCT + iPseC 
    
    # store results
    PC[[iNode]] = iPC;
  }
  return(PC)
}

# My MAP function
myMAP = function(PC){
  nodes = names(PC); n_nodes = length(nodes) # nodes
  
  MAP = list() # allocate structure
  for (i in 1:n_nodes){
    iNode = nodes[i] # node, i
    iPC = PC[[iNode]]; # table of BDe counts, i
    D_iPC = dim(iPC); nD_iPC = length(D_iPC); # dimensions of BDe table, i
    if (nD_iPC > 1){ # Nodes with parents 
      iMAP = apply(iPC, c(2:nD_iPC), function(x){ # MAP estimates
        res = x/sum(x)
        return(res)
      })
    } else { # Root nodes (nodes without parents)
      iMAP = iPC/sum(iPC) # MAP estimates
    }
    iMAP = as.table(iMAP)
    MAP[[iNode]] = iMAP
  } 
  return(MAP)
}


# One realization of all CPTs (draws from Dirichlet distribution)
myDirichletSampler = function(PC){
  require(MCMCpack) # Dirichlet implementation
  nodes = names(PC); n_nodes = length(nodes) # nodes
  pSample = list(); # allocate space for list of CPTs
  for (i in 1:n_nodes){
    iNode = nodes[i] # node, i
    iPC = PC[[iNode]]; # table of BDe counts, i
    D_iPC = dim(iPC); nD_iPC = length(D_iPC); # dimensions of BDe table, i
    if (nD_iPC > 1){ # Nodes with parents 
      i_pSample = apply(iPC, c(2:nD_iPC), function(x_i){ # sample from dirichlet
        res1 = x_i; # allocate space while preserving attributes ! (else attr will be lost)
        res0 = rdirichlet(1, alpha = x_i) # sample from dirichlet
        res1[] = res0; # updata
        return(res1)
      })
    } else { # Root nodes (nodes without parents)
      i_pSample = iPC; # allocate space while preserving attributes ! (else attr will be lost)
      res = rdirichlet(1, alpha = iPC) # sample from dirichlet
      i_pSample[] = res; # update
    }
    i_pSample = as.table(i_pSample) # store as table
    pSample[[iNode]] = i_pSample;
  }
  return(pSample)
}


# EM algorithm ------------------------------------------------------------

# E-step (unicore) - expected sufficient statistics
EM_Estep_unicore = function(FITTED, dDATA, n_iss, quiet){
  require('bnlearn')
  require('gRain')
  
  nodes = names(dDATA) # variables
  JT_base = compile(as.grain(FITTED)) # junction tree for inference (gRain)
  
  eCOUNTS <- list() # expected counts (allocation)
  Nins = nrow(dDATA) # no. of data instances
  for (n in 1:Nins){ # loop over data cases
    if (quiet == F){
      print(c('Instance', n, 'of', Nins)) # print progress
    }
    ddata_n = as.matrix(dDATA[n,]) # Data instance n
    idNA_n = which(is.na(ddata_n)) # NA index in data instance n
    idObs_n = setdiff(1:length(nodes),idNA_n) # observed elements of data instance n
    obs_n = ddata_n[idObs_n] # observations
    JT_n = setEvidence(JT_base, nodes = nodes[idObs_n], states = ddata_n[idObs_n]) # update JT based on evidence
    
    for (i in 1:length(nodes)){  # loop over variables
      nodes_ni = nodes[i] # node, i
      parents_ni = bnlearn::parents(FITTED,nodes_ni) # parents for node i
      variables_ni = c(nodes_ni, parents_ni)
      
      QG_ni = querygrain(JT_n, nodes = c(nodes_ni, parents_ni), type = "joint", exclude = FALSE) # joint distribution
      
      if (n==1){ # Update
        eCOUNTS[[nodes_ni]] = QG_ni
      } else {
        eCOUNTS[[nodes_ni]] = eCOUNTS[[nodes_ni]] + QG_ni
      }
    }
  }
  # Arrange counts as in bnlearn (the standard I use)
  eCOUNTS_bn = list()
  for (i in 1:length(nodes)){
    eCOUNTS_i = as.table(eCOUNTS[[nodes[i]]])
    nodes_i = nodes[i] # node, i
    parents_i = bnlearn::parents(FITTED,nodes_i) # parents for node i
    eCOUNTS_bn_i = aperm(eCOUNTS_i, c(c(nodes_i,parents_i)))
    eCOUNTS_bn[[nodes_i]] = eCOUNTS_bn_i
  }
  # Posterior counts
  PCprior = BDeu_prior(BN=bn.net(FITTED), dDATA=dDATA, n_iss=n_iss) # prior counts
  PC = PCprior # Initialize posterior counts to prior counts
  for (i in 1:length(nodes)){
    PC[[nodes[i]]] = PC[[nodes[i]]] + eCOUNTS_bn[[nodes[i]]] # update posterior counts by including expected counts
  }
  return(PC) # output posterior counts
}

# E-step (multicore) - expected sufficient statistics (no progress statements yet)
EM_Estep_multicore = function(FITTED, dDATA, n_iss, no_cores){
  require('bnlearn')
  require('gRain')
  require('doParallel')
  
  nodes = names(dDATA) # variables
  JT_base = compile(as.grain(FITTED)) # junction tree for inference (gRain)
  
  eCOUNTS <- list() # expected counts (allocation)
  Nins = nrow(dDATA) # no. of data instances
  
  cl <- makeCluster(no_cores, outfile="") # Define cluster
  clusterExport(cl, varlist=c('FITTED', 'dDATA'), envir=environment()) # send input to slaves
  jnk1 = clusterEvalQ(cl, c(library(gRain), library(bnlearn))) # load package to workers
  registerDoParallel(cl)
  res_list = foreach (n = 1:Nins) %dopar% { # loop over data cases
    ddata_n = as.matrix(dDATA[n,]) # Data instance n
    idNA_n = which(is.na(ddata_n)) # NA index in data instance n
    idObs_n = setdiff(1:length(nodes),idNA_n) # observed elements of data instance n
    obs_n = ddata_n[idObs_n] # observations
    JT_n = setEvidence(JT_base, nodes = nodes[idObs_n], states = ddata_n[idObs_n]) # update JT based on evidence
    
    for (i in 1:length(nodes)){  # loop over variables
      nodes_ni = nodes[i] # node, i
      parents_ni = bnlearn::parents(FITTED,nodes_ni) # parents for node i
      variables_ni = c(nodes_ni, parents_ni)
      
      QG_ni = querygrain(JT_n, nodes = c(nodes_ni, parents_ni), type = "joint", exclude = FALSE) # joint distribution
      eCOUNTS[[nodes_ni]] = QG_ni
    }
    return(eCOUNTS)
  }
  stopCluster(cl) # shut down cluster
  
  eCOUNTS = list()
  for (k in 1:length(res_list)){
    res_list_k = res_list[[k]]
    for (i in 1:length(nodes)){
      if (k==1){
        eCOUNTS[[nodes[i]]] = res_list_k[[nodes[i]]]
      } else
        eCOUNTS[[nodes[i]]] = eCOUNTS[[nodes[i]]] + res_list_k[[nodes[i]]]
    }
  }
  # Arrange counts as in bnlearn (the standard I use)
  eCOUNTS_bn = list()
  for (i in 1:length(nodes)){
    eCOUNTS_i = as.table(eCOUNTS[[nodes[i]]])
    nodes_i = nodes[i] # node, i
    parents_i = bnlearn::parents(FITTED,nodes_i) # parents for node i
    eCOUNTS_bn_i = aperm(eCOUNTS_i, c(c(nodes_i,parents_i)))
    eCOUNTS_bn[[nodes_i]] = eCOUNTS_bn_i
  }
  # Posterior counts
  PCprior = BDeu_prior(BN=bn.net(FITTED), dDATA=dDATA, n_iss=n_iss) # prior counts
  PC = PCprior # Initialize posterior counts to prior counts
  for (i in 1:length(nodes)){
    PC[[nodes[i]]] = PC[[nodes[i]]] + eCOUNTS_bn[[nodes[i]]] # update posterior counts by including expected counts
  }
  return(PC) # output posterior counts
}

# M-step - MAP parameters (See myMAP function above)

# EM algorithm (multicore)
EM_multicore = function(BN, dDATA, n_iss, no_cores, FITTED, quiet, tol, max_ite){
  require('bnlearn')
  require('gRain')
  require('doParallel')
  
  # Initialize
  PCprior=BDeu_prior(BN=BN, dDATA=dDATA, n_iss=n_iss) # prior counts
  if (is.null(FITTED)==TRUE){
    dist0 = myDirichletSampler(PCprior) # only for testing my implementation
    FITTED = custom.fit(BN, dist = dist0)
  }
  dist1 = myDirichletSampler(PCprior)
  FITTED_new = custom.fit(BN, dist = dist1)
  
  # Algorithm
  n_ite = 0
  while ((isTRUE(all.equal(FITTED,FITTED_new, tolerance = tol))==FALSE)&(n_ite < max_ite)){
    n_ite = n_ite + 1; 
    if (quiet ==FALSE){print(c('iteration', n_ite))} # progress
    if (n_ite > 1){FITTED = FITTED_new} # update
    
    PC_new = EM_Estep_multicore(FITTED=FITTED, dDATA=dDATA, n_iss=n_iss, no_cores=no_cores) # E-step
    MAP_new = myMAP(PC = PC_new) # M-step
    FITTED_new = custom.fit(BN, dist = MAP_new) # bnlearn object
  }
  return(FITTED_new)
}