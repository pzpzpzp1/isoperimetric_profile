# Medial Axis Isoperimetric Profiles

This project computes upper bounds on the isoperimetric profile using the medial axis and morphological opening.
It also contains code to compute the improved total variation lower bound. 

* Presentation of Associated Paper - [Video](https://www.youtube.com/watch?v=ynTFKEjiNcA)

### Prerequisites

Direct dependencies
* [BOOST](https://www.boost.org/) - CGAL requires this.
* [CGAL](https://www.cgal.org/download.html) - Computational Geometry Algorithms Library 
* [Multiprod](https://www.mathworks.com/matlabcentral/fileexchange/8773-multiple-matrix-multiplications-with-array-expansion-enabled) - Tensor multiplication support for matlab

IP lower bound dependencies
* [TVProfile](https://github.com/justso1/tv_profile) - Total Variation Isoperimetric Profile
* [CVX](http://cvxr.com/cvx/) - Convex problem solver for Matlab
* [Mosek](https://www.mosek.com/) - Convex problem solver. (Optional, but recommended)

## Build and Run

* Go to "ipvoronoi" folder. 
* Use cmake to generate project solution into "build" folder. (REQUIRES CGAL)
* Build "ipvoronoi" project. (this step might differ depending on OS)
* Make sure compiled executable is accessible as "isoperimetric_profile\ipvoronoi\build\Release\ipvoronoi.exe"
* Go back to "isoperimetric_profile" folder.
* Run "addpath(genpath('.'));"
* Go to "generatefigures/generateGeneralIProfiles" folder.
* Run generateGProfile.m

## Common issues

* "isotropicTotalVariation" undefined - make sure TVProfile is installed and on matlab's path. Or set params.skipTV to before calling aggregateProcessing.
* cmake isn't generating anything that can be built. - make sure CGAL and BOOST are both installed.

## Useful Functions

IP bounds
* isoperim_profile_v3 - computes IP upper bound by medial axis traversal
* morphOpenBound - computes IP upper bound by morphological opening
* getTVProfile - wrapper around TVProfile to compute IP lower bound

Visualization
* visualizeProfile - plots IP bounds of a domain
* extractAggregateIPEnvelope - combine multiple bounds into one summarizing envelope
* visualizeMedAxis - plots medial axis

Helper Functions
* getMaxInscribedCircle - gets maximum inscribed circle in a 2D domain

## Authors

* [*Paul Zhang*](https://www.csail.mit.edu/person/paul-zhang)
* [*Daryl Deford*](https://www.csail.mit.edu/person/daryl-deford)
* [*Justin Solomon*](https://people.csail.mit.edu/jsolomon/)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
