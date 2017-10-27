import numpy as np

def diff_x(f,Delta=1):
    return 0.5*(np.roll(f,-1,axis=1) - np.roll(f,1,axis=1))/Delta

def diff_y(f,Delta=1):
    return 0.5*(np.roll(f,-1,axis=0) - np.roll(f,1,axis=0))/Delta

def curl(X,Y,Dx,Dy):
    return diff_x(Y,Dx) - diff_y(X,Dy)

def norm(X,Y,Z):
    return X**2 + Y**2
