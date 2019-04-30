#!/usr/bin/env python
# coding: utf-8

# In[15]:


# Exercice 3

# Active subspace

# a) Find the actie subspace of the test function on domain [-1, 1]^10

# Import the libraries

import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C


# In[16]:


'''This exercice was based on the work of Scheidegger & Bilionis (2017).

Reference:

   Scheidegger & Bilionis (2017), 
   Machine Learning for High-Dimensional Dynamic Stochastic Economies

   https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2927400

   The example below corresponds to figure 3 from the paper.'''


# In[17]:


# The example uses the procedure for constructing the active subspace, 
# compares the Gaussian process regression on the original subspace to the GPR on the active subspace.

# Define the function to test: f(x) = exp_function(x)

# Exp(x)

def Exp_function_a(x):
    return np.exp(0.3*x[0] + 0.7*x[1])


# Calculate the gradient of Exp_function

def dExp_function_a(x):
    val = Exp_function_a(x)
    return np.column_stack([0.3*val, 0.7*val])


# In[22]:


# Evaluate the Exp_function on points in the domain = [-1,1]^2

def randOmega(N, D):
    # Use the projection matrix D
    "random points on Omega = [-1,1]^D"
    # Generate random point in the domain = [-1,1]^2
    return 2 * (np.random.rand(N, D) - 0.5)


# Define the test to be run

def test_a():
    np.random.seed(43)
    Nvals = np.array([4, 8, 16, 32])
    num_Ns = len(Nvals)
    ave_errs_as = np.inf * np.ones(num_Ns)
    ave_errs = np.inf * np.ones(num_Ns)
    
    # Construct the test points, based on the domain
    X_test = randOmega(1000, 2)
    f1_test = Exp_function_a(X_test.T)
    
    for i, N in enumerate(Nvals):
        
        # Define the training points
        X = randOmega(N, 2)
        # The Exp_function
        V = Exp_function_a(X.T)
        # Gardient of the Exp_function
        G = dExp_function_a(X.T)
        
        # Find the active subspace
        CN = (np.matmul(G.T,G)) / N
        vals, vecs = linalg.eigh(CN)
        W = vecs[:, 1]
        Y = np.matmul(X, W)
        
        # Fit Gaussian Process on active subspace
        # Using the Kernel =  C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
        kernel = RBF()
        
        # Define the Gaussian Process, based on Kernel
        gp_as = GaussianProcessRegressor(kernel)
        
        gp_as.fit(Y[:, None], V)
        m_tilde_as = gp_as.predict((np.matmul(X_test, W))[:, None])
        ave_errs_as[i] = np.sqrt(np.sum((f1_test - m_tilde_as)**2 / (f1_test**2)) / N)
        
        # Fit Gaussian Process on the input domain
        gp_as = GaussianProcessRegressor(RBF())
        gp_as.fit(X, V)
        m_tilde = gp_as.predict(X_test)
        ave_errs[i] = np.sqrt(np.sum((f1_test - m_tilde)**2 / (f1_test**2)) / N)


# Plot the figure
    fig, ax = plt.subplots(1, 2, figsize=(9, 5))
    ax[0].plot([0.9, 2.2], W, ".", ms=15)
    ax[0].set_xlabel("Input dimension")
    ax[0].set_xticks([1, 2])
    ax[0].set_ylabel("Magnitude of W")
    
    ax[1].semilogy(Nvals, ave_errs_as, ".-", label="ASGP", ms=15)
    ax[1].semilogy(Nvals, ave_errs, ".-", label="GP", ms=15)
    ax[1].legend()
    fig.tight_layout()
    plt.show()
    
    return fig


# In[23]:


# Run the example code for test_a)

test_a()


# In[ ]:


'''
b) Approximate the 10-dimensional Exp_function with GP regression
   For 10, 50, 100, 500 points randomly sampled from [-1 , 1]^10


Reference: 
    Scheidegger & Bilionis (2017), 
    Machine Learning for High-Dimensional Dynamic Stochastic Economies

    https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2927400

    The example below corresponds to figure 4 from the paper.'''


# In[29]:


# In this example, we move to a higher dimensional function test_function: ℝ^10→ℝ.


# Define the Exp_function_b function with higher dimension

def Exp_function_b(x):
    return np.exp(0.01*x[0] + 0.7*x[1] + 0.02*x[2] + 0.03*x[3] + 0.04*x[4] + 
                  0.05*x[5] + 0.06*x[6] + 0.08*x[7] + 0.09*x[8] + 0.1*x[9])

# Gardient of the Exp_function_b

def dExp_function_b(x):
    val = Exp_function_b(x)
    coefs = np.array([0.01, 0.7, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.09, 0.1])
    return val[:, None] * coefs[None, :]


# In[30]:


# Evaluate the Exp_example_b function on the domain = [-1,1]^2

# Define the domain on randomly generated points 

def randOmega(N, D):
    "random points on \Omega = [-1,1]^D"
    return 2 * (np.random.rand(N, D) - 0.5)


# Define the example function on the higher dimensional domain
def test_b():
    np.random.seed(41)
    N = 300
    X = randOmega(N, 10)
    V = Exp_function_b(X.T)
    G = dExp_function_b(X.T)
    CN = (np.matmul(G.T, G)) / N

    # Find the active subspace
    vals, vecs = linalg.eigh(CN)
    W = vecs[:, -1]
    
    # Plot the figure
    fig, ax = plt.subplots(1, 2, figsize=(9, 5))
    x = np.arange(1, 11)
    ax[0].plot(x, vals, ".", ms=15)
    ax[0].set_xlabel("Eigenvalues")
    ax[0].set_xticks(x)
    ax[0].set_ylabel("$\lambda$")
    
    ax[1].plot(x, W, ".", ms=15)
    ax[1].set_xlabel("Input dimension")
    ax[1].set_xticks(x)
    ax[1].set_ylabel("Magnitude of W")
    
    fig.tight_layout()
    plt.show()

    return fig


# In[31]:


# Run the example b)

test_b()


# In[ ]:


'''
c) Approximate the function with combination of active subspace and GP
   Compute average and max error.


Reference:
    Scheidegger & Bilionis (2017), 
    Machine Learning for High-Dimensional Dynamic Stochastic Economies

    https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2927400.'''


# In[32]:


# In this example, we plot the max absolute error on randomly generated points

# Use the old function (Exp_example_b) used in b)

def Exp_old(x):
    return np.exp(0.01*x[0] + 0.7*x[1] + 0.02*x[2] + 0.03*x[3] + 0.04*x[4] + 
                  0.05*x[5] + 0.06*x[6] + 0.08*x[7] + 0.09*x[8] + 0.1*x[9])



# Set the new function: change the domain

def Exp_function_c(x):
    return x[1]*x[2] * np.exp(0.01*x[0] + 0.7*x[1] + 0.02*x[2] + 0.03*x[3] + 0.04*x[4] + 
                  0.05*x[5] + 0.06*x[6] + 0.08*x[7] + 0.09*x[8] + 0.1*x[9])


# Define the Gardient of the function

def dExp_function_c(x):
    _Exp_old = Exp_old(x)
    val = np.atleast_1d(x[1] * x[2] * _Exp_old)
    coefs = np.array([0.01, 0.7, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.09, 0.1])
    out = val[:, None] * coefs[None, :]
    out[:, 1] += x[2] * _Exp_old
    out[:, 2] += x[1] * _Exp_old
    return out



# In[33]:


# Evaluate this function on points in the domain = [-1,1]^2

# Define the domain

def randOmega(N, D):
    "random points on \Omega = [-1,1]^D"
    
    # Randomly generated points in the domain
    return 2 * (np.random.rand(N, D) - 0.5)

# Define test_c
def test_c():
    np.random.seed(42)
    Nvals = np.array([10, 30, 100, 250, 500, 1000])
    num_Ns = len(Nvals)
    max_errs = np.inf * np.ones((num_Ns, 3))
    max_errs_gp = np.inf * np.ones(num_Ns)
    
    # Construct the test points
    X_test = randOmega(1000, 10)
    f_test = Exp_function_c(X_test.T)

    for i, N in enumerate(Nvals):
    
    # Set the training points
        X = randOmega(int(N), 10)
        V = Exp_function_c(X.T)
        G = dExp_function_c(X.T)
    
        CN = (np.matmul(G.T, G)) / N
        vals, vecs = linalg.eigh(CN)
        for d in range(3, 0, -1):

            # Find the active subspace of dimension d
            W = vecs[:, -d:]
            Y = np.matmul(X, W)

            # Fit the Gaussian Process on the active subspace
            gp_as = GaussianProcessRegressor(RBF(), n_restarts_optimizer=0)
            gp_as.fit(Y.reshape(N, d), V)
            m_tilde = gp_as.predict((np.matmul(X_test, W)).reshape(1000, d))
            max_errs[i, d-1] = np.max(np.abs(f_test - m_tilde))
        
        gp = GaussianProcessRegressor(RBF())
        gp.fit(X, V)
        m_tilde = gp.predict(X_test)
        max_errs_gp[i] = np.max(np.abs(f_test - m_tilde))
    
    fig, ax = plt.subplots(1, 2, figsize=(9, 5))
    x = np.arange(1, 11)
    ax[0].semilogy(x, vals, ".", ms=15)
    ax[0].set_xlabel("Eigenvalues")
    ax[0].set_xticks(x)
    ax[0].set_ylabel("$\lambda$")
    
    ax[1].semilogy(Nvals, max_errs, ".-", ms=15)
    ax[1].semilogy(Nvals, max_errs_gp, ".-", ms=15)
    ax[1].legend([str(i) + "d AS" for i in range(1, 4)] + ["Full GP"])
    ax[1].set_xlabel("# of points")

    fig.tight_layout()
    plt.show()
    return fig, max_errs


figure, max_errs = test_c()


# In[ ]:


'''
d) Compare the solution computed in b) and in c)


'''

