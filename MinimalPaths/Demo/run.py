#!/usr/bin/env python3

#Example runs : 
#./run.py 10. ReedsSheppForward2 1 Optimal

# System imports
import argparse
import ast
import os

# Numerics imports
from matplotlib import pyplot as plt
import matplotlib.image
import numpy as np; xp=np

# Personal imports
try: at_home = ("/Users/mirebeau/" in os.environ['CONDA_PREFIX'])
except KeyError: at_home=False
if at_home: import sys; sys.path.insert(0,"/Users/mirebeau/Dropbox/Programmes/Github/AdaptiveGridDiscretizations")

from agd import Plotting
from agd.Plotting import savefig; savefig.dirName = "."
from agd import Eikonal

if at_home: Eikonal.LibraryCall.binary_dir['FileHFM']="/Users/mirebeau/bin/HamiltonFastMarching/FileHFM/Release"
else: Eikonal.LibraryCall.binary_dir['FileHFM']="." #/Users/mirebeau/bin/HamiltonFastMarching/FileHFM/Release"

# Executable calling conventions
# with open("./workfile.txt","w") as f:
#	import platform
#	import subprocess
# 	f.write(platform.system()+"\n")
# 	if platform.system()=='Linux': dir=''
# 	elif platform.system()=='Darwin':dir='.'
# 	else:raise ValueError("Unsupported system")
# 	f.write(os.popen(os.path.join(dir,"FileHFM_ReedsSheppForward2")).read())
# #raise SystemExit(0)

# Parse input arguments
ap = argparse.ArgumentParser()
ap.add_argument("xi", type=float)
ap.add_argument("model", type=str)
ap.add_argument("nseeds", type=int)
ap.add_argument("tangent", type=str)
args = ap.parse_args()

hfmIn = Eikonal.dictIn({
	'model':args.model,
	'xi':args.xi,
	'cost':1,
})


# Input image file
im = xp.array(matplotlib.image.imread("./input_0.png"))
im = np.moveaxis(im,0,1)[:,::-1,:3] # ignore alpha channel, use cartesian coordinates

hfmIn.update({
	'walls':np.logical_not(im.sum(axis=2)>0),

	'dims':(*im.shape[:2],64), # Number of orientations as parameter ? 
	'origin':[-0.5,-0.5],
	'gridScale':1,
	'verbosity':0,
})

# Parse and normalize input points file
with open("./inpainting_data_0.txt","r") as f:
	pts = [ast.literal_eval(line.rstrip()) for line in f]
pts = np.array(pts).T
pts[1] = im.shape[1]-pts[1] # Reverse y axis
nseeds = args.nseeds
if nseeds>=pts.shape[1]//(2 if args.tangent=='Prescribed' else 1):
	print("!! Error !! Not enough points provided.")
	raise SystemExit(0)



if args.tangent=='Optimal':
	hfmIn['seeds_Unoriented'] = pts.T[:nseeds]
	hfmIn['tips_Unoriented'] = pts.T[nseeds:]
else:
	if args.tangent=='Prescribed':
		k = pts.shape[1]//2
		p = pts[:,:(2*k-1):2]
		v = pts[:,1:(2*k):2] - p
		θ = np.arctan2(v[1],v[0])
	elif args.tangent=='Random':
		p = pts
		θ = 2.*np.pi*np.random.rand(p.shape[1])
	else: 
		raise ValueError("Unexpected tangent type")

	lifted = np.concatenate((p,θ[None]),axis=0)
	hfmIn['seeds'] = lifted.T[:nseeds]
	hfmIn['tips'] = lifted.T[nseeds:]
	hfmIn['stopWhenAllAccepted']=hfmIn['tips']


# --- Output control ---
# with open("./output_workfile.txt","w") as f:
# 	f.write(f"xi={args.xi}\n")
# 	f.write(f"model={args.model}\n")
# 	f.write(f"nseeds={args.nseeds}\n")
# 	f.write(f"tangent={args.tangent}\n")
# 	f.write(f"pts={pts}.\n\n")
# 	f.write(os.popen("ls .").read())

# ------- Print the seeds and tips -------
fig = plt.figure(); plt.axis('equal'); plt.axis('off')
Plotting.imshow_ij(im)
if args.tangent=='Optimal':
	for key,color in [('seeds_Unoriented','blue'),('tips_Unoriented','red')]:
		plt.scatter(*hfmIn[key].T,color=color)
else:
	for key,color in [('seeds','blue'),('tips','red')]:
		for x,y,θ in hfmIn[key]:
			plt.arrow(x,y,20*np.cos(θ),20*np.sin(θ),color=color,head_width=10)
savefig(fig,"output_seeds_and_tips.png",dpi=100)
plt.close()


#print("Exiting before FMM run")
#raise SystemExit(0)

hfmOut = hfmIn.Run()

fig = plt.figure(); plt.axis('equal'); plt.axis('off')
Plotting.imshow_ij(im)
geodesics = (hfmOut['geodesics_Unoriented'] if args.tangent=='Optimal' 
	else hfmOut['geodesics'])
for geo in geodesics: plt.plot(*geo[:2])
savefig(fig,"output_geodesics.png",dpi=100)
plt.close()








