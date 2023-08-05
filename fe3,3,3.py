import numpy as np

def create_seq(a,b,num):
    return np.linspace(a,b,num)

print(create_seq(0.,1.,101))