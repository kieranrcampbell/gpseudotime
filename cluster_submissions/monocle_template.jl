
# coding: utf-8

# In[217]:

using Gadfly

using DataFrames

include("../bgplvm.jl"); # load in inference & plotting methods;


# In[2]:

X = readcsv("../data/X.csv")
t_gt = vec(readcsv("../data/t_gt.csv"))


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

n_iter = 500000
burn = n_iter / 2

thin = 500

# pseudotime parameters
eps = 1e-6
t =  rand(Uniform(.5 - eps, .5 + eps), n) # t_gt #  
tvar = 6.5e-3

# kernel parameters
lambda = [.5, .5]
sigma = [.5, .5]

lvar = [1e-2, 1e-2] * .9
svar = [1e-2, 1e-2] * .8

r = 1e-1 # repulsion parameter

gamma = 10.0

return_burn = true # should the burn period be returned in the traces?
cell_swap_probability = 0 # randomly swap two cells at each stage?

mh = B_GPLVM_MH(X, n_iter, burn, thin, t, tvar, lambda, lvar, 
                sigma, svar, r, return_burn, cell_swap_probability, gamma)

## write things to file
base_h5 = string("r", r)
qoi = ("tchain", "lambda_chain", "sigma_chain",
    "burn_acceptance_rate")

for q in qoi
    filename = string("../data/", base_h5, "_", q, ".csv")
    writecsv(filename, mh[q])
end

filename = string("../data/", base_h5, "_burn_thin.csv")

writecsv(filename, mh["params"]["burn_thin"])