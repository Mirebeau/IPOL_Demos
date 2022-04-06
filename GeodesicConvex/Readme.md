# Geodesic models with convexity shape prior


The IPOL demo 77777000192 [Link](https://ipolcore.ipol.im/demo/clientApp/demo.html?id=77777000192) is *inspired by* and *illustrates* the following paper:
Chen, D., Mirebeau, J.-M., Shu, M., Tai, X., and Cohen, L. D. (2021). *Geodesic Models with Convexity Shape Prior*. [arXiv Preprint](https://arxiv.org/abs/2111.00794)

More precisely:
- The IPOL demo relies on Python and C++ code, uses some of the techniques presented in the paper. (Author in charge : J.M.Mirebeau)
- The numerical experiments presented in the paper rely on Matlab code provided in the `DaChen_Matlab_Source.zip` archive. (Author in charge : Da Chen)
- **TODO** : mention the notebook


## Producing the demo files for IPOL

- The demoExtras file contains the Python code. On my machine it is generated as:
`tar -cvf demoExtras.tar -C /Users/mirebeau/Dropbox/Programmes/Github/IPOL_Demos/GeodesicConvex/Demo run.py -C /Users/mirebeau/Dropbox/Programmes/Github/AdaptiveGridDiscretizations agd`
- The HFM library is needed. For IPOL purposes, it has been compressed using the following command:
`tar -cvf hfm.tar --exclude=build --exclude=.git --exclude=*.nb --exclude=__pycache__ -C /Users/mirebeau/Dropbox/Programmes/Github HamiltonFastMarching`
-The HFM library is built on the IPOL servers, using the following part of the DDL file
        "build1":{
        "url": "https://mirebeau.github.io/IPOL_Demos/hfm.tar",
        "move":"FileHFM_ReedsSheppForward2,FileHFM_Elastica2,FileHFM_Dubins2,FileHFM_ConvexReedsSheppForward2,FileHFM_ConvexElastica2,FileHFM_ConvexDubins2",
        "construct":"cmake HamiltonFastMarching/Interfaces/FileHFM/ -DStandardModelNames:STRING='Elastica2;ReedsSheppForward2;Dubins2;ConvexElastica2;ConvexReedsSheppForward2;ConvexDubins2' -DIncludeExperimentalModels:BOOL=FALSE -DCMAKE_BUILD_TYPE:BOOL=Release; make -j4"
        }
 Related reading : 
 https://stackoverflow.com/questions/34302265/does-cmake-build-type-release-imply-dndebug
 https://stackoverflow.com/questions/12896988/passing-the-argument-to-cmake-via-command-prompt

## Testing the demo locally

In the Demo directory : `./run.py 10. Elastica2 0.05 5 2`

*Note on execution time.* About 1 minute. (The demo was intended to run on GPUs, but for now only CPUs are available on IPOL.)
