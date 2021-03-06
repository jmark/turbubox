import numpy as np
import math
import itertools

eps = np.finfo(float).eps
nit,TOL = 4,4*eps

def shiftMatrix(nrows,ncols,inmat):
    exmat = np.zeros(inmat.shape)
    for i in range(nrows,inmat.shape[0]):
        for j in range(ncols,inmat.shape[1]):
            exmat[i,j] = inmat[i-nrows,j-ncols]
    return exmat

def upAndDownScaleMat(N):
    tempMat = np.zeros((N,N))
    tempMat[0,0] = 1.0
    tempMat[1,0] = 1.0

    rfaceMat = np.zeros((2,N,N))
    for f in range(2):
        for i in range(N//2):
            rfaceMat[f] = rfaceMat[f] + shiftMatrix(2*i,f*N//2+i,tempMat)
            
    cfaceMat = np.zeros((2,N,N))
    for f in range(2):
        cfaceMat[f] = 0.5*np.transpose(rfaceMat[f])

def mk_body_centered_linspace(infimum, supremum, nNodes, withBoundaryNodes=False):
    """
    Make regular body centered linear space w/ or w/o neighboring boundary nodes.
    """

    domsize = np.abs(supremum - infimum)
    offset  = domsize / nNodes / 2

    if withBoundaryNodes:
        nNodes     = nNodes + 2
        infimum    = infimum  - offset
        supremum   = supremum + offset
    else:
        infimum    = infimum  + offset
        supremum   = supremum - offset

    return np.linspace(infimum, supremum, nNodes, endpoint=True)

def BarycentricWeight(xs,j):
    acc = 1
    for i in range(len(xs)):
        if i==j: continue
        acc *= xs[j]-xs[i]
    return 1/acc

def BarycentricWeights(xs):
    n  = len(xs)
    ws = np.empty(n)
    for j in range(n):
        ws[j] = BarycentricWeight(xs,j)     
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

def DiffMatrix(xs):
    ws = BarycentricWeights(xs)
    n = len(xs)
    M = np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            if i != j:
                M[i,j] = ws[j]/(ws[i]*(xs[i]-xs[j]))
            else:
                acc = 0
                for k in range(n):
                    if i == k: continue
                    acc += ws[k]/(xs[i]-xs[k]) 
                M[i,i] = -acc/ws[i]
    return M

def MassMatrix(ws):
    return np.diagflat(ws)

def LagrangePolynomialDerivative(xs,j,x):
    ACC = 0
    for i in range(len(xs)):
        if i == j: continue
        acc = 1
        for m in range(len(xs)):
            if m == i or m == j: continue
            acc *= (x-xs[m])/(xs[j]-xs[m])
        ACC += acc/(xs[j]-xs[i])
    return ACC

# =========================================================================== #

def LegendrePolynomialAndDerivative(N,x,doNormalize=False):
    if N == 0:
        LN, LND = 1, 0
    elif N == 1:
        LN, LND = x, 1
    else:
        LN_2,  LN_1  = 1, x
        LND_2, LND_1 = 0, 1

        for k in range(2,N+1):
            LN           = (2*k-1)/k * x * LN_1 - (k-1)/k * LN_2
            LND          = LND_2 + (2*k-1) * LN_1
            LN_2,  LN_1  = LN_1,  LN
            LND_2, LND_1 = LND_1, LND

    if doNormalize:
        LN  *= np.sqrt(N + 0.5)
        LND *= np.sqrt(N + 0.5)

    return (LN,LND)

def LegendreGaussNodesAndWeights(N):
    """ Compute the nodes (roots of the Legendre Polynomial) and weights for
        the Legendre-Gauss-Quadrature. """
    if N == 0: return (np.array([0]),np.array([2]))
    if N == 1: return (np.array([-np.sqrt(1/3),np.sqrt(1/3)]) , np.array([1,1]) )

    # nodes and weights
    xs = np.zeros(N+1)
    ws = np.zeros(N+1)
    
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
    if N == 1: return (np.array([-1,1]), np.array([1,1]))

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

def mk_mass_matrix(xs,ws):
    n = len(xs)
    M = np.empty([n,n])
    for i in range(n):
        for j in range(n):
            f = lambda x: LagrangePolynomial(xs,i,x) * LagrangePolynomial(xs,j,x)
            M[i,j] = integrate(xs,ws,f)
    return M

def mk_visual_matrix(Xs,xs):
    M = np.empty([len(Xs),len(xs)])
    for i in range(len(Xs)):
        for j in range(len(xs)):
            M[i,j] = LagrangePolynomial(xs,j,Xs[i])
    return M

def mk_corner_matrix(N):
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

# =========================================================================== #

def mk_exponent_vector(ndim,npoly):
    #return (x for x in itertools.product(range(npoly+1),repeat=ndim) if sum(x) <= npoly)
    return itertools.product(range(npoly+1),repeat=ndim)

def mk_polynome_vector(x,npoly):
    return [np.prod(np.power(x,e)) for e in mk_exponent_vector(x.shape[-1], npoly)]

def mk_polynome_matrix(xs, npoly):
    return np.array([mk_polynome_vector(x,npoly) for x in xs])

def put_row(M,nrow,row):
    tmp = M.copy()
    tmp[nrow] = row
    return tmp

def mk_interpol_matrix(xs,Xs,npoly):
    M = mk_polynome_matrix(xs,npoly)

    return np.array([
        [np.linalg.det(put_row(M,nrow,mk_polynome_vector(X,npoly)))
        for nrow in range(len(M))] for X in Xs])/np.linalg.det(M)

def mk_polynomial_interpolator(xs,Xs,npoly):
    IM = mk_interpol_matrix(xs, Xs,npoly)
    def closure(fs, domain):
        return np.dot(IM,fs)
    return closure 

# =========================================================================== #

def mk_lagrange_vector(xs,x):
    return np.array([LagrangePolynomial(xs,j,x) for j in range(len(xs))])

def mk_vandermonde_matrix(xs,ys):
    return np.array([mk_lagrange_vector(xs,x) for x in ys])

# =========================================================================== #

def mk_lagrange_interpolator_2d(xs,ys, Xs):
    def polyv(nodes,x):
        return np.array([LagrangePolynomial(nodes,j,x) for j in range(len(nodes))])

    def polyouter(x,y):
        return np.einsum('i,j->ij',polyv(xs,x),polyv(ys,y))
    
    tensors = [polyouter(*X) for X in Xs]
    
    def interpolate(fs):
        return np.array([np.sum(fs*t) for t in tensors])
        
    return interpolate

def mk_lagrange_interpolator_3d(xs,ys,zs,Xs):
    def polyv(nodes,x):
        return np.array([LagrangePolynomial(nodes,j,x) for j in range(len(nodes))])

    def polyouter(x,y,z):
        return np.einsum('i,j,k->ijk',polyv(xs,x),polyv(ys,y),polyv(zs,z))
    
    tensors = [polyouter(*X) for X in Xs]

    def interpolate(fs):
        return np.array([np.sum(fs*t) for t in tensors])
        
    return interpolate

def mk_nodes(npoly, ntype='gauss'):
    if ntype == 'gauss':
        fun = LegendreGaussNodesAndWeights
    elif ntype == 'gauss-lobatto':
        fun = LegendreGaussLobattoNodesAndWeights
    elif ntype == 'cell-centered':
        fun = lambda npoly: (mk_body_centered_linspace(-1,1,npoly+1),None)
    else:
        raise KeyError("unknown node type: '%s'" % ntype)

    nodes, _ = fun(npoly)
    return np.array(nodes)

def mk_nodes_from_to(x0,x1,npoly,ntype='gauss'):
    return x0 + (x1-x0) * (mk_nodes(npoly,ntype)+1)/2

def mk_LegendreVandermondeMatrix(xs,doNormalize=False):
    return np.array([[
                LegendrePolynomialAndDerivative(j,xs[i],doNormalize=doNormalize)[0]
            for j in range(0,len(xs))
        ] for i in range(0,len(xs))
    ])

def mk_dg2fvMat(xs,ws,M):
    """
    Args:
        xs: array of quadrature nodes
        ws: array of quadrature weights

    Returns:
        (M x len(xs)) projection matrix
    """

    xL = -0.5*np.sum(ws)
    xR =  0.5*np.sum(ws)

    # interpolation nodes within each finite volume
    Ys = [mu + xs/M for mu in mk_body_centered_linspace(xL,xR,M)]

    # projection matrix
    return np.array([np.dot(ws,mk_vandermonde_matrix(xs,ys)) for ys in Ys])

def mk_galerkinMat(xs,ws,ys,vs,M):
    """
    Args:
        xs: array of sub-element quadrature nodes
        ws: array of sub-element quadrature weights

        ys: array of main element quadrature nodes
        vs: array of main element quadrature weights

        M: number of (xs,ws)-elements

    Returns:
        tuple of M projection (len(ys) x len(xs)) matrices
    """

    xL = -0.5*np.sum(vs)
    xR =  0.5*np.sum(vs)

    # interpolation nodes within each sub-element of master element
    Ys = [mu + xs/M for mu in mk_body_centered_linspace(xL,xR,M)]

    return tuple(
        np.array([[LagrangePolynomial(ys,i,Ys[m][j])*ws[j]/vs[i]/M for j in range(len(xs))] for i in range(len(ys))])
    for m in range(M))

def mk_dg2dgMat(xs,ws,ys,vs):
    return mk_galerkinMat(xs,ws,ys,vs,M=1)
