#=
Bayesian Gaussian Process latent variable models for
pseudotime inference in single-cell RNA-seq data


kieran.campbell@sjc.ox.ac.uk=#

using Distributions
using Gadfly
using DataFrames
using StatsBase

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

function log_likelihood(X, t, lambda, sigma)
    @assert length(lambda) == length(sigma) == size(X)[2]
    n, ndim = size(X)
    ll = 0
    for i in 1:ndim
        ll += sum(logpdf(MultivariateNormal(zeros(length(t)), covariance_matrix(t, lambda[i], sigma[i])), X[:,i]))
    end
    return ll
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
            ll += log(sin(pi * abs(t[i] - t[j])))
        end
    end
    return 2 * r * ll
end

function lambda_prior(lambda, rate = 1.0)
    lp = sum(logpdf(Exponential(rate), lambda))
    # lp = 0
    return lp
end

function sigma_prior(sigma, rate = 1.0)
    # sp = sum(logpdf(InverseGamma(alpha, beta), sigma))
    sp = sum(logpdf(Exponential(rate), sigma))
    # sp = 0
    return sp
end


function acceptance_ratio(X, tp, t, lambda_prop, lambda, sigma_prop, sigma, r, s, gamma, delta)
    """ 
    Compute the acceptance ratio for 
    @param X N-by-D data array for N points in D dimensions
    @param tp Proposed pseudotime of length N
    @param t Previous pseudotime length N
    @param thetap Propose theta = [lambda, sigma]
    @param theta Previous theta
    @param r > 0 Corp parameter
    @param s Tempering parameter: (log likelihood * prior) is raised to this value
    """
    likelihood = log_likelihood(X, tp, lambda_prop, sigma_prop) - log_likelihood(X, t, lambda, sigma)
    t_prior = corp_prior(tp, r) - corp_prior(t, r)
    l_prior = lambda_prior(lambda_prop, gamma) - lambda_prior(lambda, gamma)
    s_prior = sigma_prior(sigma_prop, delta) - sigma_prior(sigma, delta)
    return s * (likelihood + t_prior + l_prior + s_prior) 
end

function couple_update_acceptance_ratio(X, t1, t2, theta1, theta2, r, s1, s2)
    h(X, t, theta, r) = log_likelihood(X, t, theta) + corp_prior(t, r) 
    return  ( (s1 - s2) * ( h(X, t2, theta2, r) - h(X, t1, theta1, r) ) )
end

function acceptance_ratio_likelihood_only(X, tp, t, thetap, theta, r, s)
    """ 
    Same as acceptance ratio except only the log_likelihood is raised to s
    
    @param s Tempering parameter: (log likelihood) is raised to this value
    """
    likelihood = log_likelihood(X, tp, thetap) - log_likelihood(X, t, theta)
    prior = corp_prior(tp, r) - corp_prior(t, r)
    #println(likelihood, " ",  prior)
    return s * likelihood + prior 
end

function couple_update_acceptance_ratio_likelihood_only(X, t1, t2, theta1, theta2, r, s1, s2)
    """ Same as couple_acceptance_ratio except only the likelihood is raised to s """
    h(X, t, theta, r) = log_likelihood(X, t, theta) 
    return  ( (s1 - s2) * ( h(X, t2, theta2, r) - h(X, t1, theta1, r) ) )
end
   
## sampling function - DO NOT USE
function sample_t(t, var)
    n = length(t)
    tprop = rand(MvNormal(t, diagm(fill(var, n))))
    return tprop
end

# TODO: propose and propose_t now same function
function propose(mean, var, lower = 0, upper = Inf)
    # sample from truncated normal of (mean, real)
    n = length(mean)
    return [rand(Truncated(Normal(mean[i], var[i]), lower, upper)) for i in 1:n]
end

function propose_t(t, var, lower = 0, upper = 1)
    #= random truncated normal for mean vector t and scalar sigma var =#
    n = length(t)
    tp = [rand(Truncated(Normal(t[i], var), lower, upper)) for i in 1:n]
    return tp
end;


function B_GPLVM_MH(X, n_iter, burn, thin, 
    t, tvar, lambda, lvar, sigma, svar, 
    r = 1, return_burn = false, cell_swap_probability = 0,
    gamma = 1.0, delta = 1.0)
    
    chain_size = int(floor(n_iter / thin)) + 1 # size of the thinned chain
    burn_thin = int(floor(burn / thin)) # size of the burn region of the thinned chain
    
    n, ndim = size(X)

    @assert ndim == 2 # for now
    @assert cell_swap_probability >= 0
    @assert cell_swap_probability <= 1
    @assert length(lambda) == length(sigma) == ndim
    @assert burn < n_iter
    @assert length(t) == n
    @assert length(lvar) == length(svar) == ndim

    ## chains
    tchain = zeros((chain_size, n))
    tchain[1,:] = t

    lambda_chain = zeros(chain_size, ndim)
    lambda_chain[1,:] = lambda

    sigma_chain = zeros(chain_size, ndim)
    sigma_chain[1,:] = sigma
    
    accepted = zeros(n_iter)

    loglik_chain = zeros(chain_size)
    prior_chain = zeros(chain_size)
    
    loglik_chain[1] = log_likelihood(X, t, lambda, sigma)
    prior_chain[1] = corp_prior(t, r)

    # alpha_chain = zeros(chain_size - 1)
    # rnd_chain = zeros(chain_size - 1)
    
    ## MH
    for i in 1:n_iter
        # proposals
        t_prop = propose_t(t, tvar)
        lambda_prop = propose(lambda, lvar)
        sigma_prop = propose(sigma, svar)

        if cell_swap_probability > 0
            if rand() < cell_swap_probability
                # swap two cells at random
                to_swap = sample(1:length(t), 2, replace = false)
                t_prop[to_swap] = t_prop[reverse(to_swap)]
            end
        end

        # calculate acceptance ratio
        alpha = acceptance_ratio(X, t_prop, t, 
                                lambda_prop, lambda, sigma_prop, 
                                sigma, r, 1, gamma, delta)

        rnd = log(rand())

        # accept - reject
        if alpha > rnd
            # accept
            accepted[i] = 1
            t = t_prop
            lambda = lambda_prop
            sigma = sigma_prop
        end
        
    
        if i % thin == 0
            # update traces
            j = int(i / thin) + 1
            tchain[j,:] = t
            lambda_chain[j,:] = lambda
            sigma_chain[j,:] = sigma
            loglik_chain[j] = log_likelihood(X, t, lambda, sigma)
            prior_chain[j] = corp_prior(t, r)
        end
    end
    
    burnt = burn_thin
    if return_burn
        burnt = 1
    end
    
    rdict = {"tchain" => tchain[burnt:end,:],
        "lambda_chain" => lambda_chain[burnt:end,:],
        "sigma_chain" => sigma_chain[burnt:end,:],
        "acceptance_rate" => (sum(accepted) / length(accepted)),
        "burn_acceptance_rate" => (sum(accepted[burnt:end]) / length(accepted[burnt:end])),
        "r" => r,
        "loglik_chain" => loglik_chain,
        "prior_chain" => prior_chain,
        "params" => {"n_iter" => n_iter,
                    "burn" => burn,
                    "thin" => thin,
                    "burn_thin" => burn_thin
            }
        }
    return rdict
end




#### Posterior predictive mean

# Posterior predictive mean given by
# $$ \mathbf{\mu_x} = K^T_* K^{-1} \mathbf{x} $$
# $$ \mathbf{\mu_y} = K^T_* K^{-1} \mathbf{y} $$
# where $K_*$ is the covariance matrix between the latent $\mathbf{t}$ and the predictive $\mathbf{t'}$, where typically $\mathbf{t'}$ is sampled as $m$ equally spaced points on the interval [0, 1].


function predict(tp, t_map, lambda_map, sigma_map, X)
    #= Returns MAP prediction of mean function given:
    @param tp Values of t at which to predict function
    @param t_map Map estimate of latent pseudotimes
    @param lambda_map Map estimate of lambda
    @param sigma_map Map estimate of sigma
    =#
    @assert length(lambda_map) == length(sigma_map) == size(X)[2]
    @assert length(t_map) == size(X)[1]
    ndim = size(X)[2]
    np = length(tp)

    Xp = fill(0.0, (np, ndim))
    for i in 1:ndim
        K_map = covariance_matrix(t_map, lambda_map[i], sigma_map[i])
        K_star_transpose = cross_covariance_matrix(tp, t_map, lambda_map[i])

        matrix_prefactor = K_star_transpose * inv(K_map)

        mu = matrix_prefactor * X[:,i]
        Xp[:,i] = vec(mu)   
    end

    return Xp
end;


#----------- Plotting functions

function plot_pseudotime_trace(mh)
    nchoose = 4
    chosen = sample(1:n, nchoose)

    df = convert(DataFrame, mh["tchain"][:, chosen])
    df[:iter] = 1:(size(df)[1])
    df = convert(DataFrame, df)
    df_melted = stack(df, [1:nchoose])
    names!(df_melted, [symbol(x) for x in ["variable", "value", "iter"]])

    return Gadfly.plot(df_melted, x = "iter", y = "value", color = "variable", Geom.line)  
end

function plot_kernel_parameter(mh, param)
    chain_name = string(param, "_chain")

    df = convert(DataFrame, mh[chain_name])
    ndim = size(df)[2]

    names!(df, [symbol(string(param, string(i))) for i in 1:ndim])
    df[:iter] = 1:(size(df)[1]) # (burn + 2)
    df_melted = stack(df, [1:2])
    return Gadfly.plot(df_melted, x = "iter", y = "value", colour = "variable", Geom.line)
end

function plot_posterior_mean(mh, tp, X)
    burn = mh["params"]["burn_thin"]
    lambda_map = mean(mh["lambda_chain"][burn:end,:], 1)
    sigma_map = mean(mh["sigma_chain"][burn:end,:], 1)

    t_map = mean(mh["tchain"][burn:end,:], 1)
    mu_p = predict(tp, t_map, lambda_map, sigma_map, X)

    return Gadfly.plot(layer(x = X[:,1], y = X[:,2], color = t_gt, Geom.point) ,
    layer(x = mu_p[:,1], y = mu_p[:,2], Geom.line(preserve_order = 1), 
    Theme(default_color=color("red"))))
end

function plot_likelihood(mh)
    df = DataFrame()
    df[:value] = mh["loglik_chain"]
    df[:iter] = 1:(size(df)[1]) # (burn + 2)
    Gadfly.plot(df, x = "iter", y = "value", Geom.line)
end

function plot_prior(mh)
    df = DataFrame()
    df[:value] = mh["prior_chain"]
    df[:iter] = 1:(size(df)[1]) # (burn + 2)
    Gadfly.plot(df, x = "iter", y = "value", Geom.line)
end


