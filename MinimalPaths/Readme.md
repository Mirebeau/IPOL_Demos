# Geodesic models with convexity shape prior


The IPOL demo 77777000191 [Link](https://ipolcore.ipol.im/demo/clientApp/demo.html?id=77777000191) is *inspired by* and *illustrates* the following paper:
Mirebeau, J.-M., Gayraud, L., Barrère, R., Chen, D., & Desquilbet, F. (2021). Massively parallel computation of globally optimal shortest paths with curvature penalization. [link](https://hal.archives-ouvertes.fr/hal-03171069/document)

However, with the caveats:
- The IPOL demo only illustrates curvature penalized with obstacles, and does not include the ground cost function found in other applications.
- The demo was intended to run on GPUs, but at the time of writing only CPUs are available on IPOL.
- Many more illustrations of curvature penalized shortest paths can be found at [this repository](https://github.com/Mirebeau/AdaptiveGridDiscretizations).


## Producing the demo files for IPOL

- The demoExtras IPOL file contains the Python code. On my machine it is generated as:
`tar -cvf demoExtras.tar -C /Users/mirebeau/Dropbox/Programmes/Github/IPOL_Demos/MinimalPaths/Demo run.py -C /Users/mirebeau/Dropbox/Programmes/Github/AdaptiveGridDiscretizations agd`
- The HFM library is needed. For IPOL purposes, it has been compressed and is built on the IPOL servers
        "build1":{
        "url": "https://mirebeau.github.io/IPOL_Demos/hfm.tar",
        "move":"FileHFM_ReedsShepp2,FileHFM_ReedsSheppForward2,FileHFM_Elastica2,FileHFM_Dubins2,",
        "construct":"cmake HamiltonFastMarching/Interfaces/FileHFM/ -DStandardModelNames:STRING='Elastica2;ReedsShepp2;ReedsSheppForward2;Dubins2' -DIncludeExperimentalModels:BOOL=FALSE -DCMAKE_BUILD_TYPE:BOOL=Release; make -j4"
        }
 Related reading : 
 https://stackoverflow.com/questions/34302265/does-cmake-build-type-release-imply-dndebug
 https://stackoverflow.com/questions/12896988/passing-the-argument-to-cmake-via-command-prompt

## Testing the demo locally (Debug only)

In the Demo directory : `./run.py 10. Elastica2 0.05 5 2`

*Note on execution time.* Possibly about 1 minute. (The demo was intended to run on GPUs, but for now only CPUs are available on IPOL.)