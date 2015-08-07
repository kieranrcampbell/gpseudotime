
# coding: utf-8

## Bayesian Gaussian Process Latent Variable models for single-cell RNA-seq pseudotime inference

##### Kieran Campbell <kieran.campbell@sjc.ox.ac.uk>

#### Setup

# In[39]:

using Distributions
using Gadfly
using PyPlot
using DataFrames


#### Functions

# In[530]:

function pairwise_distance(t)
    n = length(t)
    T = zeros((n,n))
    for i in 1:n
        for j in (i+1):n
            T[i,j] = (t[i] - t[j])^2
        end
    end
    return T + transpose(T)
end

function cross_pairwise_distance(t1, t2)
    n1 = length(t1)
    n2 = length(t2)
    T = zeros((n1, n2))
    for i in 1:n1
        for j in 1:n2
            T[i,j] = (t1[i] - t2[j])^2
        end
    end
    return T
end

function covariance_matrix(t, lambda, sigma)
    T = pairwise_distance(t)
    Sigma = exp(-lambda * T) + sigma * eye(length(t))
    return Sigma
end

function cross_covariance_matrix(t1, t2, lambda)
    T = cross_pairwise_distance(t1, t2)
    return exp(-lambda * T)
end

function log_likelihood(x, t, theta)
    lambda, sigma = theta
    return sum(logpdf(MultivariateNormal(zeros(length(t)), covariance_matrix(t, lambda, sigma)), x))
end

## electroGP

function corp_prior(t, r = 1)
    if r == 0 # uniform prior
        return 0
    end
    
    ll = 0
    n = length(t)
    for j in 1:n
        for i in (j+1):n
            ll += log(abs(sin(pi * (t[i] - t[j]))))
        end
    end
    ##println(ll, " ", r, " ", 2 * r * ll)
    return 2 * r * ll
end

function acceptance_ratio(X, tp, t, thetap, theta, r, s)
    """ 
    Compute the acceptance ratio for 
    @param X N-by-D data array for N points in D dimensions
    @param tp Proposed pseudotime of length N
    @param t Previous pseudotime length N
    @param thetap Propose theta = [lambda, sigma]
    @param theta Previous theta
    @param r > 0 Corp parameter
    @param s Tempering parameter: log likelihood is raised to this value
    """
    likelihood = log_likelihood(X, tp, thetap) - log_likelihood(X, t, theta)
    prior = corp_prior(tp, r) - corp_prior(t, r)
    #println(likelihood, " ",  prior)
    return s * (likelihood + prior) 
end

function couple_update_acceptance_ratio(X, t1, t2, theta1, theta2, r, s1, s2)
    h(X, t, theta, r) = log_likelihood(X, t, theta) + corp_prior(t, r) 
    return  ( (s1 - s2) * h(X, t2, theta2, r) + (s2 - s1) * h(X, t1, theta1, r) )
end
   
## sampling function - DO NOT USE
function sample_t(t, var)
    n = length(t)
    tprop = rand(MvNormal(t, diagm(fill(var, n))))
    return tprop
end

function propose(mean::Real, var::Real, lower = 0, upper = Inf)
    # sample from truncated normal of (mean, real)
    return rand(Truncated(Normal(mean, var), lower, upper))
end

function propose_t(t, var, lower = 0, upper = 1)
    #= random truncated normal for mean vector t and scalar sigma var =#
    n = length(t)
    tp = [rand(Truncated(Normal(t[i], var), lower, upper)) for i in 1:n]
    return tp
end;


#### Synthetic data

# In[305]:

srand(123)
lambda = 1
sigma = 1e-3
n = 100
t_gt = rand(n)
K = covariance_matrix(t_gt, lambda, sigma)
x = rand(MvNormal(zeros(n), K))
y = rand(MvNormal(zeros(n), K))
X = [x y]
Gadfly.plot(x = x, y = y, color = t_gt)


#### Metropolis-Hastings

# In[277]:

function B_GPLVM_MH(X, n_iter, burn, t, tvar, theta, theta_var, r = 1)
    n, ndim = size(X)
    @assert ndim == 2
    
    lambda, sigma = theta
    lvar, svar = theta_var
    
    tchain = zeros((n_iter + 1, n))
    tchain[1,:] = t

    theta_chain = zeros((n_iter + 1), 2)
    theta_chain[1,:] = [lambda, sigma]
    
    accepted = zeros(n_iter)

    for i in 1:n_iter
        # proposals
        t_prop = propose_t(t, tvar)
        lambda_prop = propose(lambda, lvar)
        sigma_prop = propose(sigma, svar)

        # calculate acceptance ratio
        alpha = acceptance_ratio(X, t_prop, t, 
                                [lambda_prop, sigma_prop], 
                                [lambda, sigma], r, 1)

        # accept - reject
        if alpha > log(rand())
            # accept
            accepted[i] = 1
            t = t_prop
            (lambda, sigma) = [lambda_prop, sigma_prop]
        end
        tchain[i+1,:] = t
        theta_chain[i+1,:] = [lambda, sigma]
    end
    
    rdict = {"tchain" => tchain[burn:end,:],
        "theta_chain" => theta_chain[burn:end,:],
        "acceptance_rate" => (sum(accepted) / length(accepted)),
        "burn_acceptance_rate" => (sum(accepted[burn:end]) / length(accepted[burn:end])),
        "r" => r,
        "params" => {"n_iter" => n_iter,
                    "burn" => burn}
        }
    return rdict
end;


# In[320]:

srand(120)

# sampling parameters
n_iter = 40000
burn = n_iter / 2

# pseudotime parameters
t = rand(Uniform(.499, .501), n)
tvar = .5e-3

# kernel parameters
lambda = 1
sigma = 1e-3

lvar = .5e-5
svar = .5e-10

mh = B_GPLVM_MH(X, n_iter, burn, t, tvar, [lambda, sigma], [lvar, svar], 1);

println(mh["params"])
println("Acceptance rate ", mh["acceptance_rate"])
println("Burn acceptance rate ", mh["burn_acceptance_rate"])


# In[312]:

nchoose = 20
df = convert(DataFrame, mh["tchain"][:, 1:nchoose])
df[:iter] = 1:(burn + 2)
df_melted = stack(df, [1:nchoose]);
#Gadfly.plot(df_melted, x = "iter", y = "value", color = "variable", Geom.line) ;


# In[313]:

df = convert(DataFrame, mh["theta_chain"])
names!(df, [symbol(x) for x in ["lambda", "sigma"]])
df[:iter] = 1:(burn + 2)
df_melted = stack(df, [1:2]);
Gadfly.plot(df, x = "iter", y = "lambda", Geom.line)


# In[314]:

Gadfly.plot(df, x = "iter", y = "sigma", Geom.line)


# In[315]:

Gadfly.plot(x = t_gt, y = mean(mh["tchain"], 1))


#### Posterior predictive mean

# Posterior predictive mean given by
# $$ \mathbf{\mu_x} = K^T_* K^{-1} \mathbf{x} $$
# $$ \mathbf{\mu_y} = K^T_* K^{-1} \mathbf{y} $$
# where $K_*$ is the covariance matrix between the latent $\mathbf{t}$ and the predictive $\mathbf{t'}$, where typically $\mathbf{t'}$ is sampled as $m$ equally spaced points on the interval [0, 1].

# In[318]:

m = 50
tp = linspace(0, 1, 50*m)

lambda_map = mean(mh["theta_chain"][:,1])
sigma_map = mean(mh["theta_chain"][:,2])
t_map = mean(mh["tchain"], 1)

K_map = covariance_matrix(t_map, lambda_map, sigma_map)
K_star_transpose = cross_covariance_matrix(tp, t_map, lambda_map)

matrix_prefactor = K_star_transpose * inv(K_map)

mu_x = matrix_prefactor * x
mu_y = matrix_prefactor * y;


# In[319]:

Gadfly.plot(layer(x = x, y = y, color = t_gt, Geom.point) ,
layer(x = mu_x, y = mu_y, Geom.line(preserve_order = 1), 
Theme(default_color=color("red"))))


### Parallel Tempering

# Throughout all of this we'll use $s$ for the tempering parameter since $t$ is already taken up for pseudotimes.

# In[418]:

function swap_indices(n)
    """ given n indices we can make (n-1) adjacent swaps """
    swap = sample(1:(n-1))
    return [swap, swap+1]
end;


# In[557]:

srand(123)

# temperature scales
svec = linspace(0, 1, 5) # [0, .1, .2, .5, 1] 
ns = length(svec)

r = 1

# sampling parameters
n_iter = 100000
burn = n_iter / 2

# pseudotime parameters
t = rand(Uniform(.499, .501), n)
tvar = .2e-3

# kernel parameters
lambda = 1
sigma = 1.5e-3
theta = [lambda, sigma]

lvar = .5e-5
svar = .5e-10
theta_var = [lvar, svar]

# function B_GPLVM_MH(X, n_iter, burn, t, tvar, theta, theta_var, r = 1)
n, ndim = size(X)
@assert ndim == 2

lambda, sigma = theta
lvar, svar = theta_var

"""
chains are indexed by [iteration, parameter, temperature]
"""

tchain = zeros((n_iter + 1, n, ns))


theta_chain = zeros(((n_iter + 1), 2, ns))

for i in 1:ns
    tchain[1, :, i] = t
    theta_chain[1, :, i] = theta
end

accepted = zeros(n_iter, ns)
coupled_accept = zeros(n_iter)

t = tchain[1,:,:]
lambda = vec(theta_chain[1,1,:])
sigma = vec(theta_chain[1,2,:])

t_prop = zeros(size(t)) 
lambda_prop = zeros(size(lambda))
sigma_prop = zeros(size(sigma))

alphas = zeros(ns)
saccept = fill(false, ns)

swap_position = zeros(n_iter)

for i in 1:n_iter
    if i % 10000 == 0
        println("Iter ", i)
    end
    
    ## Standard MH
    # proposals
    for j in 1:ns
        t_prop[1,:,j] = propose_t(vec(t[1,:,j]), tvar)
        lambda_prop[j] = propose(lambda[j], lvar) 
        sigma_prop[j] = propose(sigma[j], svar)  
            
        # calculate acceptance ratio
        alphas[j] = acceptance_ratio(X, vec(t_prop[1,:,j]), vec(t[1,:,j]), 
            [lambda_prop[j], sigma_prop[j]], 
        [lambda[j], sigma[j]], r, svec[j])
    end

    
    # accept - reject
    saccept = [a > log(rand()) for a in alphas]

    accepted[i,:] = 1 * saccept
    
    a_ind = find(saccept) # which temperatures accepted the update?
    t[1,:,a_ind] = t_prop[1,:,a_ind]

    lambda[a_ind] = lambda_prop[a_ind]
    sigma[a_ind] = sigma_prop[a_ind]
    
    tchain[i+1,:,:] = t
    theta_chain[i+1,:,:] = [lambda, sigma]
            
    ## Coupling update
    to_swap = swap_indices(ns) # two element array with chains to be swapped
    ## recall f & g are defined by 's'
    t1 = vec(tchain[i+1,:,to_swap[1]])
    t2 = vec(tchain[i+1,:,to_swap[2]])
    theta1 = vec(theta_chain[i+1,:,to_swap[1]])
    theta2 = vec(theta_chain[i+1,:,to_swap[2]])
    
    
    calpha = couple_update_acceptance_ratio(X, t1, t2, theta1, theta2, 
    r, svec[to_swap[1]], svec[to_swap[2]])
    # println(calpha, log(rand()))

    if(calpha > log(rand()))
        # swap chains
        coupled_accept[i] = 1
        tchain[i+1,:,to_swap[1]] = t2
        tchain[i+1,:,to_swap[2]] = t1
        theta_chain[i+1,:,to_swap[1]] = theta2
        theta_chain[i+1,:,to_swap[2]] = theta1
        
        t = tchain[i+1,:,:]
        lambda = theta_chain[i+1,1,:]
        sigma = theta_chain[i+1,2,:]
        # println("Swapped: ", string(to_swap))
        swap_position[i] = to_swap[1]
    end
end

#return rdict
# end;


# In[566]:

mh = { "tchain" => tchain[burn:end,:,end],
    "theta_chain" => theta_chain[burn:end,:,end],
        "r" => r,
    "svec" => svec,
        "params" => {"n_iter" => n_iter,
                    "burn" => burn},
    "acceptance_rate" => mean(accepted,1),
    "burn_acceptance_rate" => mean(accepted[burn:end,:], 1),
    "coupled_acceptance_rate" => mean(coupled_accept),
    "coupled_burn_acceptance_rate" => mean(coupled_accept[burn:end]),
    "swap_position" => swap_position[burn:end]
}


# In[567]:

sp = mh["swap_position"]
Gadfly.plot(x = sp[find(sp .> 0)], Geom.histogram)


# In[570]:

#Gadfly.plot(x = 1:length(mh["swap_position"]), y = mh["swap_position"], Geom.line) #;


# In[560]:

Gadfly.set_default_plot_size(14cm,10cm)
Gadfly.plot(x = t_gt, y = mean(mh["tchain"], 1))


# In[561]:

df = convert(DataFrame, mh["theta_chain"])
names!(df, [symbol(x) for x in ["lambda", "sigma"]])
df[:iter] = 1:(burn + 2)
df_melted = stack(df, [1:2]);
Gadfly.plot(df, x = "iter", y = "lambda", Geom.line)


# In[562]:

Gadfly.plot(df, x = "iter", y = "sigma", Geom.line)


# In[565]:

nchoose = 20
df = convert(DataFrame, mh["tchain"][:, 1:nchoose])
df[:iter] = 1:(burn + 2)
df_melted = stack(df, [1:nchoose]);
#Gadfly.plot(df_melted, x = "iter", y = "value", color = "variable", Geom.line)#;

