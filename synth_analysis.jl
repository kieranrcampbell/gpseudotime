

using HDF5

include("bgplvm.jl")

trace_file = "data/traces.h5";

srand(123)
lambda = 1
sigma = 1e-3
n = 100
t_gt = rand(n)
K = covariance_matrix(t_gt, lambda, sigma)
x = rand(MvNormal(zeros(n), K))
y = rand(MvNormal(zeros(n), K))
X = [x y]


# In[66]:

## Write synthetic to HDF5
h5write(trace_file, "synthetic/X", X)
h5write(trace_file, "synthetic/n", n)
h5write(trace_file, "synthetic/lambda", lambda)
h5write(trace_file, "synthetic/sigma", sigma)
h5write(trace_file, "synthetic/t_gt", t_gt)


# In[22]:

srand(120)

# sampling parameters
n_iter = 80000
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

h5write(trace_file,  "mh/tchain", mh["tchain"])
h5write(trace_file, "mh/theta_chain", mh["theta_chain"])
h5write(trace_file, "mh/n_iter",  n_iter,)

h5write(trace_file, "mh/tvar", tvar)
h5write(trace_file,  "mh/lvar", lvar)
h5write(trace_file, "mh/svar", svar)



#### Posterior predictive mean

# Posterior predictive mean given by
# $$ \mathbf{\mu_x} = K^T_* K^{-1} \mathbf{x} $$
# $$ \mathbf{\mu_y} = K^T_* K^{-1} \mathbf{y} $$
# where $K_*$ is the covariance matrix between the latent $\mathbf{t}$ and the predictive $\mathbf{t'}$, where typically $\mathbf{t'}$ is sampled as $m$ equally spaced points on the interval [0, 1].

# In[27]:

m = 50
tp = linspace(0, 1, m)

lambda_map = mean(mh["theta_chain"][:,1])
sigma_map = mean(mh["theta_chain"][:,2])
t_map = mean(mh["tchain"], 1)

K_map = covariance_matrix(t_map, lambda_map, sigma_map)
K_star_transpose = cross_covariance_matrix(tp, t_map, lambda_map)

matrix_prefactor = K_star_transpose * inv(K_map)

mu_x = matrix_prefactor * x
mu_y = matrix_prefactor * y;

h5write(trace_file, "ppm/X", [mu_x mu_y])

### Initialising with multiple chains

# In[59]:

srand(120)
nchains = 7

chain_list = Any[]

# sampling parameters
n_iter = 80000
burn = n_iter / 2

tvar = .5e-3

# kernel parameters
lambda = 1
sigma = 1e-3

lvar = .5e-5
svar = .5e-10

for i in 1:nchains
    mh = B_GPLVM_MH(X, n_iter, burn, rand(Uniform(.499, .501), n), 
    tvar, [lambda, sigma], [lvar, svar], 1)
    push!(chain_list, mh)
end



h5write(trace_file, "multi/tvar", tvar)
h5write(trace_file, "multi/lvar", lvar)
h5write(trace_file, "multi/svar", svar)
h5write(trace_file, "multi/nchains", nchains)

for i in 1:nchains
    h5write(trace_file, string("multi/chain", i,"/tchain"), chain_list[i]["tchain"])
    h5write(trace_file, string("multi/chain", i, "/theta_chain"), chain_list[i]["theta_chain"])
end

