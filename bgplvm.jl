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



function B_GPLVM_MH(X, n_iter, burn, t, tvar, theta, theta_var, r = 1, return_burn = false)
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
    burnt = burn
    if return_burn
        burnt = 1
    end
    
    rdict = {"tchain" => tchain[burnt:end,:],
        "theta_chain" => theta_chain[burnt:end,:],
        "acceptance_rate" => (sum(accepted) / length(accepted)),
        "burn_acceptance_rate" => (sum(accepted[burnt:end]) / length(accepted[burnt:end])),
        "r" => r,
        "params" => {"n_iter" => n_iter,
                    "burn" => burn}
        }
    return rdict
end

#### Posterior predictive mean

# Posterior predictive mean given by
# $$ \mathbf{\mu_x} = K^T_* K^{-1} \mathbf{x} $$
# $$ \mathbf{\mu_y} = K^T_* K^{-1} \mathbf{y} $$
# where $K_*$ is the covariance matrix between the latent $\mathbf{t}$ and the predictive $\mathbf{t'}$, where typically $\mathbf{t'}$ is sampled as $m$ equally spaced points on the interval [0, 1].


function predict(tp, t_map, lambda_map, sigma_map)
    #= Returns MAP prediction of mean function given:
    @param tp Values of t at which to predict function
    @param t_map Map estimate of latent pseudotimes
    @param lambda_map Map estimate of lambda
    @param sigma_map Map estimate of sigma
    =#
    K_map = covariance_matrix(t_map, lambda_map, sigma_map)
    K_star_transpose = cross_covariance_matrix(tp, t_map, lambda_map)

    matrix_prefactor = K_star_transpose * inv(K_map)

    mu_x = matrix_prefactor * x
    mu_y = matrix_prefactor * y

    return [mu_x mu_y]
end