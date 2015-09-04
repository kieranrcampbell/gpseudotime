
# coding: utf-8

# In[217]:

using Gadfly
using HDF5
using DataFrames

include("../bgplvm.jl"); # load in inference & plotting methods;


# In[2]:

X = h5read("../data/embeddings.h5", "monocle/embedding")
t_gt = h5read("../data/embeddings.h5", "monocle/pseudotime");


# In[3]:



# In[4]:

# remove cells less than 0 on x
to_keep = X[:,1] .> 0
X = X[to_keep,:]
t_gt = t_gt[to_keep]

x = X[:,1] ; y = X[:,2]
xs = (x - mean(x)) / sqrt(var(x))
ys = (y - mean(y)) / sqrt(var(y))

# In[52]:


#plot(x = X[:,1], y = X[:,2], colour = t_gt)
#X_old = X
X = [xs ys];
#size(X)


#### Using large kernel parameters

# In[234]:

srand(123)

n = size(X)[1]

n_iter = 200
burn = n_iter / 2

thin = 1

# pseudotime parameters
eps = 1e-6
t =  rand(Uniform(.5 - eps, .5 + eps), n) # t_gt #  
tvar = 6.5e-3

# kernel parameters
lambda = [.5, .5]
sigma = [.5, .5]

lvar = [1e-2, 1e-2] * .9
svar = [1e-2, 1e-2] * .8

r = 1e-3 # repulsion parameter

gamma = 1.0

return_burn = true # should the burn period be returned in the traces?
cell_swap_probability = 0 # randomly swap two cells at each stage?

mh = B_GPLVM_MH(X, n_iter, burn, thin, t, tvar, lambda, lvar, 
                sigma, svar, r, return_burn, cell_swap_probability, gamma)

## write things to file
base_h5 = string("gamma", gamma)
qoi = ("tchain", "lambda_chain", "sigma_chain",
    "burn_acceptance_rate")
for q in qoi
    h5write("../data/gamma_tests.h5",
        string(base_h5, "/", q), mh[q])
end

h5write("../data/gamma_tests.h5",
    string(base_h5, "/burn_thin"), mh["params"]["burn_thin"])

