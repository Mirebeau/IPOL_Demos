#!/usr/bin/env python3

# Example input
#./run.py 10. Elastica2 0.05 5 2

# System imports
import argparse
import ast
import os

# Numerics imports
from matplotlib import pyplot as plt
import matplotlib.image
import numpy as np; xp=np
from scipy import ndimage
npeigh = np.linalg.eigh

# Personal imports
try: at_home = ("/Users/mirebeau/" in os.environ['CONDA_PREFIX'])
except KeyError: at_home=False
if at_home: import sys; sys.path.insert(0,"/Users/mirebeau/Dropbox/Programmes/Github/AdaptiveGridDiscretizations")

from agd import Plotting
from agd.Plotting import savefig; savefig.dirName = "."
from agd import Eikonal
from agd import LinearParallel as lp
from agd import Metrics

if at_home: Eikonal.LibraryCall.binary_dir['FileHFM']="/Users/mirebeau/bin/HamiltonFastMarching/FileHFM/Release"
else: Eikonal.LibraryCall.binary_dir['FileHFM']="." #/Users/mirebeau/bin/HamiltonFastMarching/FileHFM/Release"

show=False # show some figures along execution, debug mode

# Parse input arguments
ap = argparse.ArgumentParser()
ap.add_argument("xi", type=float)
ap.add_argument("model", type=str)
ap.add_argument("_lambda", type=float)
ap.add_argument("rho", type=float)
ap.add_argument("sigma", type=float)
args = ap.parse_args()

dx = 1/(2*args.rho)

# Input image file
im = xp.array(matplotlib.image.imread("./input_0.png"))
im = np.moveaxis(im,0,1)[:,::-1,:3] # ignore alpha channel, use cartesian coordinates

# --- Input verification control ---

# with open("./output_workfile.txt","w") as f:
#     f.write("Hi there.\n")
#     f.write(f"xi={args.xi}\n")
#     f.write(f"model={args.model}\n")
#     f.write(f"lambda={args._lambda}\n")
#     f.write(f"rho={args.rho}\n")
#     f.write(f"sigma={args.sigma}.\n\n")
#     f.write(os.popen("ls .").read())

# ------- Eikonal solver preliminary setup ----------

hfmIn = Eikonal.dictIn({
    'model':args.model,
    'xi':args.xi,
    'cost':1,
    'dims':(*im.shape[:2],64),
    'gridScale':1.,
    'origin':(-0.5,-0.5),
})


# ----------- Structure tensor -----------
def array_map(f,x,fdim):
    """
    Input :
    - f function, which applies to arrays with fdim dimensions
    - x array with more than fdim dimensions
    Output : 
    - f mapped to array, seen as a collection of arrays with fdim dimensions
    """
    xdim = x.ndim-fdim; assert xdim>=0
    xshape,fshape = x.shape[:xdim],x.shape[xdim:]
    out = xp.array([f(xi) for xi in x.reshape((-1,*fshape))])
    out = np.reshape( out, xshape+out[0].shape)
    return out

def gaussian_filter(x,*args,fdim=2,**kwargs): # Applies to multi-channel images
    return array_map(lambda xi:ndimage.gaussian_filter(xi,*args,**kwargs),x,fdim)

def StructureTensor(u,σ=2.,ρ=5.,dx=1,fdim=2):
    """
    Computes the structure tensor of u,
    with noise scale σ and feature scale ρ.
    Applies to fdim dimensional arrays.
    - dx : grid scale.
    """
    # Compute grad uσ
    eye = np.eye(fdim).astype(int)
    duσ = [gaussian_filter(u,σ,order=e,fdim=fdim)/dx for e in eye]
    
    # Self outer product and averaging 
    S = gaussian_filter(lp.outer_self(duσ),ρ,fdim=fdim)
    
    # sum over channels
    if fdim<u.ndim: S = np.sum(S.reshape((fdim,fdim,-1,*u.shape[-fdim:])),axis=2)
    return S


# ----------- Metric construction -----------

def eigh(S): 
    """Similar to np.linalg.eigh, but with geometry first"""
#    μ,v = xp.linalg.eigh(np.moveaxis(S,(0,1),(-2,-1)))
     # cupy.linalg.eigh buggy (10.3) -> Transfering back and forth to cpu
    μ,v = [xp.array(e) for e in npeigh(np.moveaxis(S,(0,1),(-2,-1)))]
    return np.moveaxis(μ,-1,0), np.moveaxis(v,(-1,-2),(0,1))


def EdgeMetric(μ,α=0.01,λ=1.,κ=2.):
    """
    A metric enhancing motion along edges, inspired by 
    Weickert's edge enhancing diffusion filter,
    based on the structure tensor. 
    Input :
     - μ eigenvalues of the structure tensor (μ[0] <= μ[1])
     - α lower bound on anisotropy
     - λ sensitivity to small details
     - κ sensitivity to anisotropy
    Output : eigenvalues of the metric
    """
    ν = np.maximum(1e-10,μ[1] - κ*μ[0]) / λ # Edge detector. Non-negative. 
    λ0 = np.maximum(α,1.-np.exp(-3.314/ν**4)) # Small diffusion if a feature is detected
    λ1 = np.ones_like(ν) # Unit diffusion in general
    return λ0,λ1 


#-----------------------Cost function--------------------------

S = StructureTensor(np.moveaxis(im,-1,0),args.sigma,args.rho,dx=dx)
riemann = Metrics.Riemann.from_mapped_eigenvalues(S,lambda μ: EdgeMetric(μ,λ=args._lambda))


_,_,aθ = hfmIn.Axes()
cost = np.moveaxis(xp.array([riemann.norm([np.cos(θ),np.sin(θ)]) for θ in aθ]),0,-1)
hfmIn['cost'] = cost

if show:
    plt.figure()
    X = xp.array(np.meshgrid(*hfmIn.Axes()[:2],indexing='ij'))
    plt.contourf(*X,np.min(cost,axis=2))
    plt.colorbar()
    plt.axis('equal')
    plt.show()
    plt.close()

# -------------- boundary points --------------
#Find the tangent direction to the contour, 
with open("./inpainting_data_0.txt","r") as f:
    boundary_pts = [ast.literal_eval(line.rstrip()) for line in f]
boundary_pts = np.array(boundary_pts).T
boundary_pts[1] = im.shape[1]-boundary_pts[1] # Reverse y axis

bpts = np.round(boundary_pts).astype(int)
_,tangents = eigh(S[...,bpts[0],bpts[1]])
tangents=tangents[:,0] # Keep only the eigenvector associated with the small eigenvalue

# Find the clockwise orientation
assert len(boundary_pts)>=2
interior_pt = np.mean(boundary_pts,axis=1)
outward_dir = boundary_pts - interior_pt[...,None]
tangents *= np.sign( lp.det([outward_dir,tangents]))

# Construct lifted point
angles = np.mod(np.arctan2(tangents[1],tangents[0]),2*np.pi)
lifted_pts = np.concatenate((boundary_pts,angles[None]),axis=0)


if show:
    plt.figure()
    Plotting.imshow_ij(im,cmap='Greys')
    for p,v in zip(boundary_pts.T,tangents.T): plt.arrow(*p,*10*v,color='r',head_width=10)
    plt.scatter(*boundary_pts,color='r')
    plt.scatter(*interior_pt,color='b')
    plt.show()
    plt.close()

# ------------------ obstacles ---------------
def angular_wall(hfmIn):
    sθ,tθ=hfmIn['seed'][2],hfmIn['tip'][2]
    if tθ>sθ:tθ-=2*np.pi
    mθ = ((sθ+tθ)/2)%(2*np.pi)
    
    dθ = hfmIn.gridScales[2]
    _,_,θs = hfmIn.Grid()
    return np.abs(θs-mθ)<=dθ

# --------- Setup and run the FMM --------
iseed = 0
itip = (iseed+1)%lifted_pts.shape[1]
hfmIn['seed'] = lifted_pts[:,iseed]
hfmIn['tip'] = lifted_pts[:,itip]
hfmIn['walls'] = angular_wall(hfmIn)
hfmIn['stopWhenAllAccepted']=hfmIn['tips']

#print("Exiting before FMM run")
#raise SystemExit(0)

hfmOut = hfmIn.Run()

fig = plt.figure(); plt.axis('off')
Plotting.imshow_ij(im,cmap='Greys')
for p,v in zip(boundary_pts.T,tangents.T): plt.arrow(*p,*10*v,color='b',head_width=10)
plt.scatter(*boundary_pts,color='b')
for geo in hfmOut['geodesics']: plt.plot(*geo[:2],color='r')
savefig(fig,"output_geodesics.png",dpi=100)

if show:
    plt.show()
    plt.close()
