
import numpy as np
import matplotlib.pyplot as plt

i1 = [2,1,100]
i2 = [5,3,80]
i3 = [0, 4, 150]
i4 = [1,0,0]
i5 = [0,1,0]

rr = np.arange(0, i1[2]/i1[0], 1)

def y(x): 
    return ((i1[2] - i1[0]*x)/i1[1])
def y2(x):
    return ((i2[2] - i2[0]*x)/i2[1])
def y3(x): 
    return ((i3[2] - i3[0]*x)/i3[1])


def optim(i1,i2,i3,i4,i5):
    plt.plot(rr, y(rr), rr, y2(rr), rr, y3(rr))
    plt.xlim(0,50)
    plt.ylim(0,50)
    plt.title("From optimization_a")
    

    


