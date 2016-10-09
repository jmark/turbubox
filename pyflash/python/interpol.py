import numpy as np
import math

eps = np.finfo(float).eps
nit,TOL = 4,4*eps

def Weigth(xs,j):
    acc = 1
    for i in range(len(xs)):
        if i==j: continue
        acc *= xs[j]-xs[i]
    return 1/acc

def Weights(xs):
    n  = len(xs)
    ws = np.empty(n)
    for j in range(n):
        ws[j] = Weigth(xs,j)     
    return ws

def LagrangePolynomial(xs,j,x):
    acc = 1
    for i in range(len(xs)):
        if i==j: continue
        acc *= (x-xs[i])/(xs[j]-xs[i])
    return acc

def BarycentricPolynomial(xs,ws,fs,x):
    numerator,denominator = 0,0
    for j in range(len(xs)):
        diff = x-xs[j]
        if abs(diff) <= eps: return fs[j]
        numerator   += fs[j]*ws[j]/diff
        denominator +=       ws[j]/diff
    return numerator/denominator

def DiffMatrix(xs,ws):
    n = len(xs)
    M = np.empty([n,n])
    for i in range(n):
        for j in range(n):
            if i != j:
                M[i,j] = ws[j]/ws[i]/(xs[i]-xs[j])
            else:
                acc = 0
                for k in range(n):
                    if i == k: continue
                    acc += ws[k]/(xs[i]-xs[k]) 
                M[i,i] = -acc/ws[i]
    return M

def MassMatrix(ws):
    return np.diagflat(ws)

# =========================================================================== #

def LegendrePolynomialAndDerivative(N,x):
    if N == 0: return (1,0)
    if N == 1: return (x,1)

    LN_2,LN_1,LN    = 1,x,None
    LND_2,LND_1,LND = 0,1,None

    for k in range(2,N+1):
        LN    = (2*k-1)/k * x * LN_1 - (k-1)/k * LN_2
        LND   = LND_2 + (2*k-1) * LN_1
        LN_2,LN_1   = LN_1,LN
        LND_2,LND_1 = LND_1,LND

    return (LN,LND)

def LegendreGaussNodesAndWeights(N):
    """ Compute the nodes (roots of the Legendre Polynomial) and weights for
        the Legendre-Gauss-Quadrature. """
    if N == 0: return (0,2)
    if N == 1: return ( [-np.sqrt(1/3),np.sqrt(1/3)] , [1,1] )

    # list of nodes and weights
    xs,ws = np.empty(N+1),np.empty(N+1)
    
    for j in range(math.floor((N+1)/2)):
        # make initial guess for the jth's node
        xs[j] = -math.cos((2*j+1)/(2*N+2) * math.pi)     
        for k in range(0,nit):
            # Newton's method for finding the root
            LN1,LND1 = LegendrePolynomialAndDerivative(N+1,xs[j])
            Delta    = -LN1/LND1
            xs[j]   += Delta
            if abs(Delta) <= TOL * abs(xs[j]): break

        # get final optimal values for Legendre Polynomial
        LN1,LND1 = LegendrePolynomialAndDerivative(N+1,xs[j])
        ws[j]    = 2/(1-xs[j]**2)/LND1**2
        # utilize symmetry
        xs[N-j]  = -xs[j]
        ws[N-j]  = ws[j]
        
    # consider the middle point if there is one (always zero)
    if N % 2 == 0:
        LN1,LND1 = LegendrePolynomialAndDerivative(N+1,0.0)
        xs[N//2] = 0
        ws[N//2] = 2/LND1**2

    return (xs,ws)

# =========================================================================== #

def qAndLEvaluation(N,x):
    LN_2,LN_1,LND_2,LND_1  = 1,x,0,1

    for k in range(2,N+1):
        LN    = (2*k-1)/k * x * LN_1 - (k-1)/k * LN_2
        LND   = LND_2 + (2*k-1) * LN_1
        LN_2,LN_1   = LN_1,LN
        LND_2,LND_1 = LND_1,LND

    k    = N+1   
    LN1  = (2*k-1)/k * x * LN - (k-1)/k * LN_2
    LND1 = LND_2 + (2*k-1) * LN_1
    q    = LN1  - LN_2
    qD   = LND1 - LND_2

    return (q,qD,LN)

def LegendreGaussLobattoNodesAndWeights(N):
    if N == 1: return ( [-1,1], [1,1] )

    xs,ws = np.empty(N+1),np.empty(N+1)

    xs[0],xs[N] = -1,1
    ws[0],ws[N] = (2/N/(N+1),) * 2

    for j in range(math.floor((N+1)/2)):
        xs[j] = - math.cos((j+1/4)*math.pi/N - 3/8/N/math.pi/(j+1/4))
        
        for k in range(nit+1):
            q,qD,LN = qAndLEvaluation(N,xs[j])
            Delta   = -q/qD
            xs[j]  += Delta
            if abs(Delta) <= TOL*abs(xs[j]): break

        q,qD,LN = qAndLEvaluation(N,xs[j])
        xs[N-j] = -xs[j]
        ws[N-j] = ws[j] = 2/N/(N+1)/LN**2

    if N % 2 == 0:
        q,qD,LN  = qAndLEvaluation(N,0.0)
        xs[N//2] = 0
        ws[N//2] = 2/N/(N+1)/LN**2

    return (xs,ws)

# =========================================================================== #

def integrate(xs,ws,f):
    acc = 0
    for j in range(len(xs)):
        acc += ws[j] * f(xs[j]) 
    return acc

# =========================================================================== #

def mkMassMatrix(xs,ws):
    n = len(xs)
    M = np.empty([n,n])
    for i in range(n):
        for j in range(n):
            f = lambda x: LagrangePolynomial(xs,i,x) * LagrangePolynomial(xs,j,x)
            M[i,j] = integrate(xs,ws,f)
    return M

def mkVisualMatrix(Xs,xs):
    M = np.empty([len(Xs),len(xs)])
    for i in range(len(Xs)):
        for j in range(len(xs)):
            M[i,j] = LagrangePolynomial(xs,j,Xs[i])
    return M

def mkEdgeMatrix(N):
    B = np.zeros([N,N])
    B[0,0] = -1
    B[-1,-1] = 1

    return B

# =========================================================================== #

def f1(N=None):
    return lambda x: np.cos(x)

def f2(N=None):
    return lambda x: 1/(1+x**2)

def f3(N=None):
    return lambda x: x**(2*N-2)

def f4(N=None):
    return lambda x: x**(2*N)

def f5(N=None):
    return lambda x: x**(2*N+2)
