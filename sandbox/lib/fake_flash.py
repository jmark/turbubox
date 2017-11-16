# This class is a mess and only for testing, but it works. No docs here. RTFSC.
class FakeFile:
    def __init__(self, domain, boxdims, fillby='constant', blockdims=None):
       	self.domain     = np.array(domain)
        self.gridsize   = np.array(boxdims).astype(np.int)
        self.grid       = np.array([[0,0,0], self.gridsize-1])
        self.domainsize = np.abs(self.domain[1]-self.domain[0])
        self.cellsize   = self.domainsize / self.gridsize

        if blockdims:
            self.blocksize = np.array(blockdims)
        else:
            self.blocksize = None

        self.fillby = fillby
 
    @staticmethod
    def plateau(x,y,z):
        if np.abs(x) <= 0.2 and np.abs(y) <= 0.2 and np.abs(z) <= 0.2:
            return 1
        return 0

    @staticmethod
    def wiggle(X,Y,Z):
        return np.sin(4 * 2*np.pi * X) + np.sin(5 * 2*np.pi * Y) + np.sin(6 * 2*np.pi * Z)

    def data(self, dname):
        dom = self.domain.transpose()
        grd = self.gridsize
        ret = np.zeros(grd)
        fillby = self.fillby

        def exp3D(A,sigma,p,x,y,z):
            return A * np.exp(-((x-p[0])**2 + (y-p[1])**2 + (z-p[2])**2)/sigma/2.)

        X,Y,Z = np.meshgrid(*tuple(ulz.mk_body_centered_linspace(d[0],d[1],s) for d,s in zip(dom, grd)), indexing='ij')

        if dname == 'dens':

            D1 =  1.0
            D2 =  2.0

            ret = 0.1 + 0.9 * np.exp(-((Y-0.5)**2)/2/0.01)
            #ret = D1 * np.ones_like(X)
            #ret[np.where((6/16 <= Y) * (Y <= 10/16))] = D2 
            #ret[np.where((6/16 <= Y) * (Y <= 10/16))] = D2 

            # --------------------------------------------------------------- #
            # if fillby == 'constant':
            #     ret[:] = 1.0

            # elif fillby == 'planeX':
            #     ret[:] = X

            # elif fillby == 'planeX+':
            #     ret = np.where(
            #         (2/8-1/grd[0]<X) * (X<6/8+1/grd[0]) * \
            #         (2/8-1/grd[1]<Y) * (Y<6/8+1/grd[1]) * \
            #         (2/8-1/grd[2]<Z) * (Z<6/8+1/grd[2]), X,0*X)

            # elif fillby == 'planeXYZ':
            #     ret = np.where(
            #         (1/8-1/grd[0]<X) * (X<7/8+1/grd[0]) * \
            #         (1/8-1/grd[1]<Y) * (Y<7/8+1/grd[1]) * \
            #         (1/8-1/grd[2]<Z) * (Z<7/8+1/grd[2]), X+Y+Z,0*X)

            # elif fillby == 'plane+wiggle':
            #     ret = np.where(
            #         (1/16-1/grd[0]<X) * (X<15/16+1/grd[0]) * \
            #         (1/16-1/grd[1]<Y) * (Y<15/16+1/grd[1]) * \
            #         (1/16-1/grd[2]<Z) * (Z<15/16+1/grd[2]),
            #         X+Y+Z + 0.5*self.wiggle(X+1/16,Y+1/16,Z+1/16),0*X)

            # elif fillby == 'gaussianXYZ':
            #     ret = np.exp(-((X-0.5)**2 + (Y-0.5)**2 + (Z-0.5)**2)/2/0.02)

            # elif fillby == 'stepsXYZ':
            #     ret = X + 100*X*(np.sin(4*2*np.pi*X) + 0.2*np.cos(20*2*np.pi*X) + 0.1*np.sin(20*2*np.pi*X))**2 

            # else:
            #     raise NotImplementedError('unknow fillby: %s' % fillby)

        elif dname == 'pres':
            # constant
            #ret[:] = 1.0
            # plane
            #ret = X+Y

            # --------------------------------------------------------------- #
            # shear flow
            ret = np.ones_like(X)
            #ret[np.where((0.375 <= Y) * (Y <= 0.625) * (0.375 <= Z) * (Z <= 0.625))] = 1.0
            #ret[np.where((6/16 <= Y) * (Y <= 10/16))] = 5.0
            #ret[np.where((Y < 5/16) + (11/16 < Y))] = -5.0

            # gaussian2d
            # ret = 1 + 1 * 10**1 * np.exp(-((X-0.5)**2 + (Y-0.5)**2)/2/0.02)

            # gaussian3d
            #return 1 + 1 * 10**1 * np.exp(-(X**2 + Y**2 + Z**2)/2/0.02)

            # sine3d
            #return np.abs(np.sin(10*(X**2 + Y**2 + Z**2)))

            # cube2d
            # return 1 + 2 * np.vectorize(plateau)(X,Y,0)

            # wiggle2d
            #ret = 1 + 10 * np.abs(np.sin((X-Y)/2/np.pi * 200) + np.sin((X+Y)/2/np.pi * 100))
        
            # sin2d
            #ret = 200 * (np.sin(X/np.pi/2) + np.cos(Y/np.pi/2))

            #ret = np.sin(2*np.pi*X)
   
            # sin2d
            #return X+Y+Z

        elif dname == 'velx':
            ## smashing balls
            #V0 = -0.8
            #V1 = -2 * V0 
            #s1 = 0.01
            #p1 = [0.5,0.5,0.5] 

            #ret = V0 * np.ones_like(X) + exp3D(V1, s1, p1, X,Y,Z)

            # --------------------------------------------------------------- #
            ## beam
            #ret = 0.1 * (2 * np.random.rand(*X.shape) - 1)
            #ret[np.where((0.375 <= Y) * (Y <= 0.625) * (0.375 <= Z) * (Z <= 0.625))] = 1.0
            #ret[np.where((Y < 23/64) + (41/64 < Y))] = -5.0

            U1 = -0.2
            U2 =  0.4

            #ret = U1 * np.ones_like(X)
            #ret[np.where((6/16 <= Y) * (Y <= 10/16))] = U2 

            ret = U1 * np.ones_like(X) + U2 * np.exp(-((Y-0.5)**2)/2/0.01)

        elif dname == 'vely':
            U1 =  0.01
            ret = U1 * np.sin(4*np.pi*X)

        elif dname in 'velz magx magy magz'.split():
            ret[:] = 0.0

        else:
            raise KeyError('%s not found!' % dname) 

        return ret


