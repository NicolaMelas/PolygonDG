# PolygonDg Library

**Author :** Nicola Melas

**Mailto :** nicola.melas@mail.polimi.it

## The project

#### PolygonDg
PolygonDg is a library for the approximation of Elastodynamics problems with Discontinuous
Galerkin finite element methods on polygonal grids.  
It makes use of a Polygenerator to build the mesh. The generator is originally built in MATLAB. The code supports also
the use of Triangles, thanks to the flexible use of pointers.
It is written in C++14.  
Developed with g++ 7.4.0, tested with g++ 7.1.0, 5.5.0, clang++ 6.0 and clang++ 4.0.

#### Examples
The folder `main` contains examples of use of the library and some verifications
about something done inside the code. In particular there are scripts that export
the mesh to be viewed in paraview, to check the correctness of the generators.
Then there are two scripts, one for polygons, one for triangles, to check
the calculation of the convergence order of the method, with respect to the
refinement of the mesh.

## How to install

#### Requirments
* The library **Eigen** is used for linear algebra, it does not have to be compiled
  since it is a header only library.  
  It is required the **version 3.3.3** or higher. It is suggested the version 3.3.7.  
  See http://eigen.tuxfamily.org/index.php?title=Main_Page for more information.

* The library **CGAL** is used for triangulation of a polygon, it needs to be
  linked.  
  It is suggested the version 4.9 or higher.  
  See https://www.cgal.org/index.html for more information.

* The library **boost** is used for some operations on polygons such as intersections
  of polygons, rather than to compute an ordering algorithm.  
  It is suggested the version 1.72.0 or higher.  
  See https://www.boost.org/ for more information.

* The library **muparserx** is used for evaluating functions written inside
  strings. It's useful also if you combine it with GetPot to avoid recompilation.
  It is suggested the version 4.08 or higher.  
  See https://beltoforion.de/article.php?a=muparserx for more information.

* The library **qhull** is used to construct a voronoi tesselletion of a domain.
  This library is used by MATLAB also.
  It is suggested the version 2019.1. **libqhullcpp** and **libqhull_r** are used.  
  See http://www.qhull.org/ for more information.

* The library **GetPot** is used to parse the tests commaind line arguments, it is not
  necessary for building the library but it is needed if you want to run the tests
  or the examples provided about the convergence analisys. It is only one header
  file, so you do not need to compile anything.  
  See http://getpot.sourceforge.net/ for more information.

* The tool **Doxygen** is used to generate the documentation.  
  It is suggested the version 1.8.11 or higher.  
  See http://www.stack.nl/~dimitri/doxygen/index.html for more information.

#### Installation of the library
First of all you have to edit the `Makefile.inc` in this folder and insert all the
variables that are needed, in particular where the libraries header files are and
where the libraries are, if not in the default folders, and the directory
 where you want to install the library (`PolygonDg_PATH`).

Then enter the folder `libPolygonDg` and type:
```shell
make all RELEASE=yes
```
This builds both a static and a dynamic version of the library, builds the executables
of the tests and the documentation.

If you want to copy it in another folder type:
```shell
make install
```
you can install the library in the directory specified with `PolygonDg_PATH`.

The option `RELEASE=yes` is recommended in order to enable all the optimizations
and get the full performance from the code. Without it you compile in the debug mode.

Typing `make help` you can get some information about other kinds of commands.

If the library has been successfully built, you should find in the folder `libPolygonDg/lib`
the two libraries, static and dynamic.  
In the folder `libPolygonDg/doc` there should be the documentation, in HTML and LaTeX format.
To compile the latex sources with pdflatex go into the folder `libPolygonDg/doc/latex` and type `make`.  
In the main page of the documentation you can find a tutorial for the use of the library.  
In the folder `libPolygonDg/bin` there should be
the executables of the tests, that run with the configuration files `namefile.txt`,
in which you can set the parameters and sources of the problem.
If you have build the dynamic version of the library before building the tests (for
example if you have used `make all`), then the executables have been linked with the
dynamic one, so remember to tell the loader where the libraries are editing
the environmental variable `LD_LIBRARY_PATH`.

If you want to uninstall the library removing all the files that had been copied
in `PolygonDg_PATH`, type:
```shell
make uninstall
```

#### Examples
All the examples are already built if you have used `make all`, otherwise, after
the creation of the library type:
```shell
make test RELEASE=yes
```
in order to compile the examples. Again remember to specify `RELEASE=yes` for the optimization.

Then you can run the executables and you can modify the .txt files to change parameters
of the problem and some other options.
