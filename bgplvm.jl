#=
Bayesian Gaussian Process latent variable models for
pseudotime inference in single-cell RNA-seq data


kieran.campbell@sjc.ox.ac.uk=#

using Distributions


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
    @param s Tempering parameter: (log likelihood * prior) is raised to this value
    """
    likelihood = log_likelihood(X, tp, thetap) - log_likelihood(X, t, theta)
    prior = corp_prior(tp, r) - corp_prior(t, r)
    #println(likelihood, " ",  prior)
    return s * (likelihood + prior) 
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


function B_GPLVM_MH(X, n_iter, burn, thin, 
    t, tvar, theta, theta_var, r = 1, return_burn = false)
    
    chain_size = int(floor(n_iter / thin)) + 1 # size of the thinned chain
    burn_thin = int(floor(burn / thin)) # size of the burn region of the thinned chain
    
    n, ndim = size(X)
    @assert ndim == 2
    
    lambda, sigma = theta
    lvar, svar = theta_var
    
    ## chains
    tchain = zeros((chain_size, n))
    tchain[1,:] = t

    theta_chain = zeros(chain_size, 2)
    theta_chain[1,:] = [lambda, sigma]
    
    accepted = zeros(n_iter)

    loglik_chain = zeros(chain_size)
    prior_chain = zeros(chain_size)
    
    loglik_chain[1] = log_likelihood(X, t, [lambda, sigma])
    prior_chain[1] = corp_prior(t)
    
    ## MH
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
        
    
        if i % thin == 0
            # update traces
            j = int(i / thin) + 1
            tchain[j,:] = t
            theta_chain[j,:] = [lambda, sigma]
            loglik_chain[j] = log_likelihood(X, t, [lambda, sigma])
            prior_chain[j] = corp_prior(t)
        end
    end
    
    burnt = burn_thin
    if return_burn
        burnt = 1
    end
    
    rdict = {"tchain" => tchain[burnt:end,:],
        "theta_chain" => theta_chain[burnt:end,:],
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
    K_map = covariance_matrix(t_map, lambda_map, sigma_map)
    K_star_transpose = cross_covariance_matrix(tp, t_map, lambda_map)

    matrix_prefactor = K_star_transpose * inv(K_map)

    x = X[:,1]
    y = X[:,2]

    mu_x = matrix_prefactor * x
    mu_y = matrix_prefactor * y

    return [mu_x mu_y]
end;


#----------- Plotting functions

function plot_pseudotime_trace(mh)
    nchoose = 50
    chosen = sample(1:n, nchoose)

    df = convert(DataFrame, mh["tchain"][:, chosen])
    df[:iter] = 1:(size(df)[1])
    df = convert(DataFrame, df)
    df_melted = stack(df, [1:nchoose])
    names!(df_melted, [symbol(x) for x in ["variable", "value", "iter"]])

    return Gadfly.plot(df_melted, x = "iter", y = "value", color = "variable", Geom.line)  
end

function plot_kernel_parameter(mh, param)
    df = convert(DataFrame, mh["theta_chain"])

    names!(df, [symbol(x) for x in ["lambda", "sigma"]])
    df[:iter] = 1:(size(df)[1]) # (burn + 2)
    df_melted = stack(df, [1:2]);
    return Gadfly.plot(df, x = "iter", y = param, Geom.line)
end

function plot_posterior_mean(mh, tp, X)
    burn = mh["params"]["burn_thin"]
    (lambda_map, sigma_map) = mean(mh["theta_chain"][burn:end,:], 1)
    t_map = mean(mh["tchain"][burn:end,:], 1)
    mu_p = predict(tp, t_map, lambda_map, sigma_map, X)

    return Gadfly.plot(layer(x = X[:,1], y = X[:,2], color = t_gt, Geom.point) ,
    layer(x = mu_p[:,1], y = mu_p[:,2], Geom.line(preserve_order = 1), 
    Theme(default_color=color("red"))))
end

function plot_likelihood(mh)
    df = convert(DataFrame, mh["loglik_chain"])
    df[:iter] = 1:(size(df)[1]) # (burn + 2)
    return Gadfly.plot(df, x = "iter", y = param, Geom.line)
end

function plot_prior(mh)
    df = convert(DataFrame, mh["prior_chain"])
    df[:iter] = 1:(size(df)[1]) # (burn + 2)
    return Gadfly.plot(df, x = "iter", y = param, Geom.line)
end