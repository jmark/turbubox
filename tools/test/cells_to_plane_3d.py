#!/usr/bin/env python

import sys
import numpy as np
import titanic
import interpolate as itpl

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.family': 'Monospace', 'font.size': 8})
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from mpl_toolkits.mplot3d.art3d import Line3DCollection
import ulz

def vnorm(v):
    return np.sqrt(np.sum(v*v))

fp = sys.argv[1]
fh = titanic.File(fp)

levels = fh.get('/mesh/levels')
morton = fh.get('/mesh/morton')

cells = np.transpose(fh.get('/data/states/dens'),(0,3,2,1))

#p = np.array([0.5,0.0,0.0],dtype=np.double)
#u = np.array([0.0,1.0,0.0],dtype=np.double)
#v = np.array([0.0,0.0,1.0],dtype=np.double)

#p = np.array([0.0,0.5,0.0],dtype=np.double)
#u = np.array([1.0,0.0,0.0],dtype=np.double)
#v = np.array([0.0,0.0,1.0],dtype=np.double)

#p = np.array([0.0,0.0,0.5],dtype=np.double)
#u = np.array([1.0,0.0,0.0],dtype=np.double)
#v = np.array([0.0,1.0,0.0],dtype=np.double)

#p = np.array([0.0,0.0,0.0],dtype=np.double)
#u = np.array([1.0,0.0,0.0],dtype=np.double)
#v = np.array([0.0,1.0,1.0],dtype=np.double)

#p = np.array([0.0,0.0,0.1],dtype=np.double)
#u = np.array([1.0,0.0,0.0],dtype=np.double)
#v = np.array([0.0,1.0,0.1],dtype=np.double)

#p = np.array([0.0,0.0,0.0],dtype=np.double)
#u = np.array([1.0,0.0,0.0],dtype=np.double)
#v = np.array([0.0,1.0,0.5],dtype=np.double)

#p = np.array([0.0,0.0,0.0],dtype=np.double)
#u = np.array([1.0,0.0,0.0],dtype=np.double)
#v = np.array([0.0,1.0,0.1],dtype=np.double)

#S = 1./np.sqrt(2)
#Q = 0.5*(1-S)
#p = np.array([0.0, Q , Q ],dtype=np.double)
#u = np.array([1.0,0.0,0.0],dtype=np.double)
#v = np.array([0.0, S , S ],dtype=np.double)

#p = np.array([0.0, 0.0, 0.0],dtype=np.double)
#u = np.array([1.0, 0.0, 0.0],dtype=np.double)
#v = np.array([0.0, 1.0, 1.0],dtype=np.double)

#p = np.array([0.0, 0.0, 0.0],dtype=np.double)
#u = np.array([1.0, 0.0, 0.5],dtype=np.double)
#v = np.array([0.0, 1.0, 0.5],dtype=np.double)

def interpol_3d(p,u,v,cells,shape=2*(2*128,)):
    image = np.zeros(shape,dtype=np.double)
    itpl.cells_to_plane_3d(levels,morton,cells,image,p,u,v,method='linear')
    return image

def get_grid3d(p,u,v):
    edgecount = itpl.plane_morton_to_coords(levels,morton,p,u,v)
    edges = np.zeros((edgecount,2,3))
    edgecount = itpl.plane_morton_to_coords(levels,morton,p,u,v,edges)
    return edges

dpi = 150

#extent = [0,1,0,1]
#extent = [0,1,0,1.01]
#extent = [-0.5,1.5,-0.5,1.5]
#extent = [0,np.linalg.norm(u),0,np.linalg.norm(v)]

#segments = edges[:,:,1:3]
#segments = edges[:,:,0:3:2]
#segments = edges[:,:,0:2]

#plt.gca().add_collection(matplotlib.collections.LineCollection(segments,color='white',linewidths=0.5,linestyles='solid'))

# plt.imshow(
#     data,
#     extent=extent,
#     vmin = 0.0, vmax = 3.0,
#     cmap='cubehelix',
#     #cmap='viridis',
#     interpolation=None,
#     origin='lower left'
# )
# plt.colorbar()

def mk_facecolors(data,cmap='cubehelix',cmin=None,cmax=None):
    if cmin is None: cmin = np.min(data)
    if cmax is None: cmax = np.max(data)

    norm = matplotlib.colors.Normalize(cmin, cmax)
    m = plt.cm.ScalarMappable(norm=norm,cmap=cmap)
    m.set_array([])
    return m.to_rgba(data)

fcmin = 0.0
fcmax = 3.0

fig = plt.figure(figsize=(720/dpi, 720/dpi), dpi=dpi)

ax = fig.gca(projection='3d')
#ax.view_init(elev=10.0,azim=-90.0)
#ax.view_init(elev=10.0,azim=-45.0)
ax.view_init(elev=10.0,azim=-120.0)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)
ax.set_zlim(0.0, 1.0)

# -------------------------------------------------------------------------- #

if 1:
    p = np.array([ 0.0,  0.0,  0.0],dtype=np.double)
    u = np.array([ 1.0,  0.0,  0.5],dtype=np.double)
    v = np.array([ 0.0,  1.0,  0.5],dtype=np.double)

    data = interpol_3d(p,u,v,cells)
    edges = get_grid3d(p,u,v)

    s = ulz.mk_body_centered_linspace(0,1,data.shape[0])
    t = ulz.mk_body_centered_linspace(0,1,data.shape[1])

    ss,tt = np.meshgrid(s,t,indexing='ij')

    # create vertices for a rotated mesh (3D rotation matrix)
    X = p[0] + ss*u[0] + tt*v[0]
    Y = p[1] + ss*u[1] + tt*v[1]
    Z = p[2] + ss*u[2] + tt*v[2]

    fc = mk_facecolors(data,cmin=fcmin,cmax=fcmax)

    ax.plot_surface(X,Y,Z, rstride=4, cstride=4, facecolors=fc, vmin=fcmin, vmax=fcmax, shade=False)
    ax.add_collection(Line3DCollection(edges,color='black',linewidths=0.05,linestyles='solid'))

# -------------------------------------------------------------------------- #

if 0:
    p = np.array([ 0.0,  0.0,  0.0],dtype=np.double)
    u = np.array([ 0.0,  0.0,  1.0],dtype=np.double)
    v = np.array([ 1.0,  1.0,  0.0],dtype=np.double)

    data = interpol_3d(p,u,v,cells)
    edges = get_grid3d(p,u,v)

    s = ulz.mk_body_centered_linspace(0,1,data.shape[0])
    t = ulz.mk_body_centered_linspace(0,1,data.shape[1])

    ss,tt = np.meshgrid(s,t,indexing='ij')

    # create vertices for a rotated mesh (3D rotation matrix)
    X = p[0] + ss*u[0] + tt*v[0]
    Y = p[1] + ss*u[1] + tt*v[1]
    Z = p[2] + ss*u[2] + tt*v[2]

    fc = mk_facecolors(data,cmin=fcmin,cmax=fcmax)

    ax.plot_surface(X,Y,Z, rstride=4, cstride=4, facecolors=fc, vmin=fcmin, vmax=fcmax, shade=False)
    ax.add_collection(Line3DCollection(edges,color='black',linewidths=0.02,linestyles='solid'))

# -------------------------------------------------------------------------- #

if 0:
    p = np.array([ 0.5,  0.5,  0.0],dtype=np.double)
    u = np.array([ 0.0,  0.0,  1.0],dtype=np.double)
    v = np.array([ 0.5, -0.5,  0.0],dtype=np.double)

    data = interpol_3d(p,u,v,cells)
    edges = get_grid3d(p,u,v)

    s = ulz.mk_body_centered_linspace(0,1,data.shape[0])
    t = ulz.mk_body_centered_linspace(0,1,data.shape[1])

    ss,tt = np.meshgrid(s,t,indexing='ij')

    # create vertices for a rotated mesh (3D rotation matrix)
    X = p[0] + ss*u[0] + tt*v[0]
    Y = p[1] + ss*u[1] + tt*v[1]
    Z = p[2] + ss*u[2] + tt*v[2]

    fc = mk_facecolors(data,cmin=fcmin,cmax=fcmax)

    ax.plot_surface(X,Y,Z, rstride=4, cstride=4, facecolors=fc, vmin=fcmin, vmax=fcmax, shade=False)
    ax.add_collection(Line3DCollection(edges,color='black',linewidths=0.1,linestyles='solid'))

# -------------------------------------------------------------------------- #

if 0:
    p = np.array([ 0.5,  0.5,  0.0],dtype=np.double)
    u = np.array([ 0.0,  0.0,  1.0],dtype=np.double)
    v = np.array([ 0.5, -0.5,  0.0],dtype=np.double)

    data = interpol_3d(p,u,v,cells)
    edges = get_grid3d(p,u,v)

    s = ulz.mk_body_centered_linspace(0,1,data.shape[0])
    t = ulz.mk_body_centered_linspace(0,1,data.shape[1])

    ss,tt = np.meshgrid(s,t,indexing='ij')

    # create vertices for a rotated mesh (3D rotation matrix)
    X = p[0] + ss*u[0] + tt*v[0]
    Y = p[1] + ss*u[1] + tt*v[1]
    Z = p[2] + ss*u[2] + tt*v[2]

    fc = mk_facecolors(data,cmin=fcmin,cmax=fcmax)

    ax.plot_surface(X,Y,Z, rstride=4, cstride=4, facecolors=fc, vmin=fcmin, vmax=fcmax, shade=False)
    ax.add_collection(Line3DCollection(edges,color='black',linewidths=0.1,linestyles='solid'))

# -------------------------------------------------------------------------- #

plt.title('density blob: cake cut')

fig.tight_layout()
plt.savefig('test.png', bbox_inches='tight', dpi=dpi)
plt.close()
