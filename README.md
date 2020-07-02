# Medial Axis Isoperimetric Profiles

This project computes upper bounds on the isoperimetric profile using the medial axis.
It also contains code to compute the improved total variation lower bound. 

### Prerequisites

Direct dependencies
* [BOOST](https://www.boost.org/) - CGAL requires this.
* [CGAL](https://www.cgal.org/download.html) - Computational Geometry Algorithms Library 
* [Multiprod](https://www.mathworks.com/matlabcentral/fileexchange/8773-multiple-matrix-multiplications-with-array-expansion-enabled) - Tensor multiplication support for matlab

IP lower bound dependencies
* [TVProfile](https://github.com/justso1/tv_profile) - Total Variation Isoperimetric Profile
* [CVX](http://cvxr.com/cvx/) - Convex problem solver for Matlab
* [Mosek](https://www.mosek.com/) - Convex problem solver. (Optional, but recommended)

## Build

* Go to ipvoronoi folder. 
* Use cmake to generate project solution. (REQUIRES CGAL)
* Build ipvoronoi project. (this step might differ depending on OS)
* Make sure compiled executable is accessible as "isoperimetric_profile\ipvoronoi\build\Release\ipvoronoi.exe"
* Go back to isoperimetric_profile folder.
* Run addpath(genpath('.'));
* Go to generatefigures/generateGeneralIProfiles folder.
* Run generateGProfile.m

## Authors

* [*Paul Zhang*](https://www.csail.mit.edu/person/paul-zhang)
* [*Daryl Deford*](https://www.csail.mit.edu/person/daryl-deford)
* [*Justin Solomon*](https://people.csail.mit.edu/jsolomon/)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
