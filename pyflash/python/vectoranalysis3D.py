import numpy as np

def diff_x(f,Delta=1):
    return (np.roll(f,-1,axis=0) - np.roll(f,1,axis=0))/2./Delta

def diff_y(f,Delta=1):
    return (np.roll(f,-1,axis=1) - np.roll(f,1,axis=1))/2./Delta

def diff_z(f,Delta=1):
    return (np.roll(f,-1,axis=2) - np.roll(f,1,axis=2))/2./Delta

def curl(X,Y,Z,Dx,Dy,Dz):
    dX = (diff_y(Z,Dy) - diff_z(Y,Dz))
    dY = (diff_z(X,Dz) - diff_x(Z,Dx))
    dZ = (diff_x(Y,Dx) - diff_y(X,Dy))
    
    return (dX,dY,dZ)

def norm(X,Y,Z):
    return X**2 + Y**2 + Z**2
