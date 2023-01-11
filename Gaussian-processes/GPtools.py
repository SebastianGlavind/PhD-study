import numpy as np

##################################################################
# KERNEL FURNCTIONS
##################################################################

# Kernel function: Squared exponential (SE)
def se_kern(x1, x2, sig2_f, l2_f):
    '''
    x - input array [NxM]
    y - output vector [Nx1]
    sig2_f, l2_f - kernel hyper-parameters [scalars]
    '''
    sqdist = np.sum(x1**2,1).reshape(-1,1) + np.sum(x2**2,1) - 2*np.dot(x1, x2.T) # reshape(-1,1): from one-dim. to two-dim. array.
    K = sig2_f * np.exp( -0.5 * sqdist / l2_f )
    return(K)
def Ky_kern(x1, x2, sig2_y, sig2_f, l2_f):
    '''
    x - input array [NxM]
    y - output vector [Nx1]
    sig2_y - noise variance
    sig2_f, l2_f - kernel hyper-parameters [scalars]
    '''
    K = se_kern(x1, x2, sig2_f, l2_f)
    Ky = K + np.eye(len(x1))*sig2_y
    return(Ky)

##################################################################
# GP MARGINAL LIKELIHOOD
##################################################################

# Marginal log-likelihood 
def logMargLike(x, y, sig2_y, sig2_f, l2_f):
    '''
    GP marginal likelihood, see Rasmussen & Williams (2006), Algo. 2.1
    x - input array [NxM]
    y - output vector [Nx1]
    sig2_y - noise variance
    sig2_f, l2_f - kernel hyper-parameters [scalars]
    '''
    N = x.shape[0]
    
    # Apply the kernel function to our training points
    Ky = Ky_kern(x, x, sig2_y, sig2_f, l2_f)
    Ly = np.linalg.cholesky(Ky)

    # inverse of covariance
    invLy = np.linalg.inv(Ly)
    invKy = invLy.T.dot(invLy) # invK = invL^T invL

    # alpha
    alpha = invKy.dot(y)
    
    LML = - 0.5 * y.T.dot( alpha ) - np.sum( np.log( np.diag(Ly) ) ) - 0.5 * N * np.log(2*np.pi)
    return( LML )

##################################################################
# GP POSTERIOR INFERENCE
##################################################################

# Compute posterior density at test points
def compute_pos(Xtr, ytr, Xte, sig2_y, sig2_f, l2_f):
    '''
    GP posterior density, see Rasmussen & Williams (2006), Algo. 2.1
    Xtr - training inputs [NxM]
    Xte - test inputs [NxM]
    ytr - training outputs [Nx1]
    sig2_y - noise variance
    sig2_f, l2_f - kernel hyper-parameters [scalars]
    '''
    nX_tr = Xtr.shape[0]
    nX_te = Xte.shape[0]

    # Evaluate kernel
    K = se_kern(Xtr, Xtr, sig2_f, l2_f)
    K_s = se_kern(Xtr, Xte, sig2_f, l2_f)
    K_ss = se_kern(Xte, Xte, sig2_f, l2_f)
    
    # Apply Cholesky decomposition
    Ly = np.linalg.cholesky(K + np.eye(nX_tr)*sig2_y)

    # inverse of covariance
    invLy = np.linalg.inv(Ly)
    invKy = invLy.T.dot(invLy) # invK = invL^T invL

    # alpha
    alp_y = invKy.dot(ytr)

    # predictive mean
    mu_pos = np.dot(K_s.T, alp_y)

    # predictive variance
    vv_m = np.linalg.solve(Ly,K_s)
    cov_pos = K_ss - vv_m.T.dot(vv_m)
        
    return(mu_pos, cov_pos)

# Sample posterior at test points
def sample_pos(mu_pos, cov_pos, n_samp):
    eps = 1e-8
    Nte = mu_pos.shape[0]
    # Draw samples from the posterior at our test points.
    L_pos = np.linalg.cholesky(cov_pos + eps*np.eye(Nte))
    f_pos = mu_pos + np.dot(L_pos, np.random.normal(size=(Nte, n_samp)))
    
    return(f_pos)
    
    