import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

"""
Scheidegger & Bilionis (2017), 
Machine Learning for High-Dimensional Dynamic Stochastic Economies
https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2927400
This example is slightly different than what S&B present in that I plot the maximum 
absolute error on the 1,000 random testing points instead of 
what they refer to as the average error.
"""

#Adjusting the dimension of the functions and coeficient to 8 dimensions
def test_old(x):
    return np.exp(0.01 * x[0] + 0.7 * x[1] + 0.02 * x[2] + 0.03 * x[3] + 0.04 * x[4] +
                  0.05 * x[5] + 0.06 * x[6] + 0.08 * x[7] + 0.09)


def test_function(x):
    return x[1]* x[2] *np.exp(0.01 * x[0] + 0.7 * x[1] + 0.02 * x[2] + 0.03 * x[3] + 0.04 * x[4] +
                                0.05 * x[5] + 0.06 * x[6] + 0.08 * x[7])


def dtest_function(x):
    _test_old = test_old(x)
    val = np.atleast_1d(_test_old)
    coefs = np.array([0.01, 0.7, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08])
    out = val[:, None] * coefs[None, :]
    out[:, 1] += x[2] * _test_old
    out[:, 2] += x[1] * _test_old
    return out


# Evaluating this function on points in \Omega = [-1,1]^8
def randOmega(N, D=8):
    "random points on \Omega = [-1,1]^D"
    return 2 * (np.random.rand(N, D) - 0.5)


def test_case():
    np.random.seed(42)
    print("Length of Nvals = number of iternations, with 100 sample points added per step. i.e. setting up the function for 10 iterations of 100 sample points per iteration")
    Nvals = np.array([100, 100, 100, 100, 100, 100, 100, 100, 100, 100])
    num_Ns = len(Nvals)
    print("Length of Nvals becomes also the number of iteration steps (= 10) used in the model")
    max_errs = np.inf * np.ones((num_Ns, 3))
    max_errs_gp = np.inf * np.ones(num_Ns)
    avg_errs_gp = np.inf * np.ones(num_Ns)

    print("construct test points, setting N = 100 test points and D = 8 dimensions")
    X_test = randOmega(100, 8)
    f_test = test_function(X_test.T)

    for i, N in enumerate(Nvals):
        # training points
        X = randOmega(int(N), 8)
        V = test_function(X.T)
        G = dtest_function(X.T)

        CN = (G.T @ G) / N
        vals, vecs = linalg.eigh(CN)
        for d in range(3, 0, -1):
            # find active subspace of dimension d
            W = vecs[:, -d:]
            Y = X @ W

            # fit GP on active subspace
            gp_as = GaussianProcessRegressor(RBF(), n_restarts_optimizer=8)
            gp_as.fit(Y.reshape(N, d), V)
            m_tilde = gp_as.predict((X_test @ W).reshape(100, d))
            max_errs[i, d - 1] = np.max(np.abs(f_test - m_tilde))

        gp = GaussianProcessRegressor(RBF())
        gp.fit(X, V)
        m_tilde = gp.predict(X_test)
        # Calculating the max error for the GP
        max_errs_gp[i] = np.max(np.abs(f_test - m_tilde))
        #Calculating the same for avg error according to the formula in part b), and using these values for part b)
        avg_errs_gp[i] = (1/N) * np.sum(np.abs(f_test - m_tilde))

#Adjusting the plot to fit 8 dimensions
    fig, ax = plt.subplots(1, 2, figsize=(9, 5))
    x = np.arange(1, 9)
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

    print("b) Plotting the average and maximum errors according to the formulae in part b:")
    plt.semilogy(Nvals, max_errs_gp, ".-", ms = 15)
    plt.title("Maximum errors")
    plt.xlabel("# of points")
    plt.show()

    plt.semilogy(Nvals, avg_errs_gp, ".-", ms=15)
    plt.title("Average errors")
    plt.xlabel("# of points")
    plt.show()
    print("As expected we get a vertical line. This line shows us the difference in errors generated by the model when we run the model 10 times over the same amount of sample points")

    print("is the gradient vector, hence, we plot and put out the value for the gradient at the 10th iteration step")
#
    plt.plot(x, G[10], ".", ms = 15)
    plt.xlabel("Dimension")
    plt.ylabel("Partial derivative")
    plt.show()
    print("The gradient vector for the 10th generation: ")
    print(G[10])

    print("For this exercise we utilise the partial derivative at the 10th iteration, where G is the matrix of partial derivatives\
    Then we construct a matrix C whose eigenvalues can help us decided the active subspaces in the 10th iteration. We see\
    From our plot that the active subspace is present at the 2nd and 3rd dimension for the 10th iteration.")


    return fig, max_errs

figure, max_errs = test_case()

# d)
def dtest_example(x):
    val = test_function(x)
    coefs = np.array([0.01, 0.7, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08])
    return val[:, None] * coefs[None, :]

def example():
    np.random.seed(42)
    # Setting N = 100 to give our model 100 sample points
    N = 100
    # Setting randOmega(N, D) where the number of dimensions is equal to D = 8
    X = randOmega(N, 8)
    G = dtest_example(X.T)
    CN = (G.T @ G) / N

    # find active subspace
    vals, vecs = linalg.eigh(CN)
    W = vecs[:, -1]

    fig, ax = plt.subplots(1, 2, figsize=(9, 5))
    # Sorting the eigenvalues
    x = np.arange(1, 9)
    ax[0].plot(x, vals, ".", ms=15)
    ax[0].set_xlabel("Sorted Eigenvalues")
    ax[0].set_xticks(x)
    ax[0].set_ylabel("$\lambda$")

    ax[1].plot(x, W, ".", ms=15)
    ax[1].set_xlabel("Input dimension")
    ax[1].set_xticks(x)
    ax[1].set_ylabel("Magnitude of W")

    fig.tight_layout()
    plt.show()

    return fig


example();
print("Answer to question 4 d): We can see from the sorted eigenvalues that there is only 1\
      active subspace in this model. And from the input dimensions we see that the active\
      subspace we observe is in the 2nd dimension")