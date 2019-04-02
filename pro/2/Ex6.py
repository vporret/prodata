
#%reset -f

#Part 1 : 

def f(x):
    return 0.25*x**4


class Central:
    def __init__(self, f, h=0.001):
        self.f, self.h = f, h
        
    def __call__(self, x):
        f, h = self.f, self.h
        return (f(x + h) - f(x - h)) / 2/h
    
# make function-like object df
df = Central(f)    # df(x) computes the derivative of f(x) approximately:
for x in (1, 5, 10):
    df_value = df(x) # approx value of derivative of f at point x
    exact = x**3
    # exact value of derivative
    print("f’(%d)=%g (error=%.2E)" % (x, df_value, exact-df_value))
    
# Part 2 : 

from math import log

print('\n%-15s %6s' % ('', 'x=10')) #create headline

for h in (0.5, 0.1, 1E-3, 1E-5, 1E-7, 1E-9, 1E-11):  # for every h in the problem set
    df = Central(log, h)  # to use log(x+h) in the class (log instead of f)
    
    print('h =', h)
    
    print('%-11s | %10.2E' % ("f' (exact)", 1 / 10))    #derivative of ln(x) = 1/2*(x)
    
    print('%-11s | %10.2E' % ("f' (approx)", df(10)))   # take 10 as in problem set

    print('%-11s | %10.2E'  % ("Error : ", 1 / 10 - df(10))) # difference between the two


