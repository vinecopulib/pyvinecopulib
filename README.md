pyvinecopulib
------------------

[![Build Status](https://travis-ci.com/MattClarkson/pyvinecopulib.svg?branch=master)](https://travis-ci.com/MattClarkson/pyvinecopulib)
[![Build Status](https://ci.appveyor.com/api/projects/status/5pm89ej732c1ekf0/branch/master)](https://ci.appveyor.com/project/MattClarkson/cmakecatchtemplate)


Purpose
-------

This is a demo project to demonstrate a reasonable structure for [CMake](https://cmake.org/) based projects,
that use [CTest](https://cmake.org/) to run unit tests via [Catch](https://github.com/catchorg/Catch2).


Credits
-------

This project was developed as a teaching aid for UCL's ["Research Computing with C++"](http://rits.github-pages.ucl.ac.uk/research-computing-with-cpp/)
course developed by [Dr. James Hetherington](http://www.ucl.ac.uk/research-it-services/people/james)
and [Dr. Matt Clarkson](https://iris.ucl.ac.uk/iris/browse/profile?upi=MJCLA42).


Main Features
-------------

The main idea of this project is to try and help a C++ algorithm developer write their algorithm in a
C++ library and then deploy the C++ library in a variety of scenarios. 

The main features provided are:

 1. A Meta-Build, also known as a SuperBuild, to optionally download and build any of the following: Boost, Eigen, FLANN, OpenCV, glog, gflags, VTK, PCL and [ArrayFire](https://arrayfire.com/). All of these can be left OFF and skipped. This results in a top-level build folder containing the compiled dependencies, and then a sub-folder containing the compiled code of this project.
 2. A single library into which you can provide your main algorithms.
 3. Unit tests, using Catch, and run with CTest, so you can ensure correctness and enable regression testing of your functionality.
 4. A single command line application, to give the end user a minimalist runnable program.
 5. KWStyle config, so you can check for consistent code style, when you have KWStyle installed on your system.
 6. CppCheck config, so you can check for some performance, style and correctness issues, when you have CppCheck installed on your system.
 7. Doxygen config, so you can generate documentation via a ```make docs``` or DOCS task in Visual Studio.
 8. If your code is open-source, you can register with a Continuous Integration service, so this project provides Travis and Appveyor examples.
 9. Basic examples of how to create a Qt+VTK, Qt+OpenGL or QML+VTK user interface, ensuring the VTK render engine works in Qt or QML framework, on Windows, Linux and Mac.
10. CPack setup to produce installers for the GUI apps.
11. An example of the CMake required to build python interfaces to your C++ code, using ```boost::python```.
12. An example of the CMake required to build python interfaces to your C++ code, using [pybind11](https://github.com/pybind/pybind11), with credit to [this example](https://github.com/pybind/cmake_example).
13. An example of the CMake required to export a C-style module into [Unity](https://unity3d.com/).
14. Support for OpenMP, which is passed through to FLANN and OpenCV.
15. Support for CUDA, which is passed through to FLANN, OpenCV and PCL.
16. Support for MPI, which by default sets up the C++ libraries.
17. If doing Boost.Python and OpenCV, an example of passing a numpy ndarray to OpenCV, computing something, and returning a cv::Mat as a numpy ndarray, thanks to [Gregory Kramida's pyboostcvconverter](https://github.com/Algomorph/pyboostcvconverter).
18. Support for Python Wheels, thanks to [Matthew Brett's multibuild](https://github.com/matthew-brett/multibuild).


Use Cases
---------

The feature list above is quite diverse, and in reality most developers 
won't want all of the above features. So, while you could attempt all of these, 
within the same project, we envisage the following project types:

* A C++ library in C++ command line interfaces. 
* A C++ library in a C++ user interface, using Qt, or QML.
* A C++ library in a Python module, or perhaps a Python GUI, where the Python is developed separately, outside of this project.
* A C++ library in an environment such as Unity, developed outside of this project.

So, the aim for the C++ developer would be to write a small algorithm library, and get the code out there
into the hands of their users. This project can accommodate all of the above, and there are CMake
options to switch these options OFF/ON.

Furthermore, once you have a successful build for your use-case, you could 
delete the bits of code that you don't need.  We could have made a different
repository for each of the use-cases, but then there would be a lot of code
duplication and overlap. So, for now, its all one repository. If you want to
streamline your project, you could:

* Take a look in the ```Code``` sub-folder. Remove the directories you do not need,
and change ```Code/CMakeLists.txt``` accordingly.

* Search the top level CMakeLists.txt for code that looks like
```mpAddSomething``` and ```mpIncludeSomething``` and easily chop out 3rd party libraries
you do not need. 

We believe its easier to remove code you don't need than it is to
build up a lot of CMake code, based on hours of searching the internet.
At least the CMake code here has been tested ... to some degree.


Basic Build Instructions
------------------------

This project itself can be built if you just want to evaluate it. In Linux terms that
would be:
``` cmake
git clone --recursive https://github.com/MattClarkson/pyvinecopulib
mkdir pyvinecopulib-Build
cd pyvinecopulib-Build
ccmake ../pyvinecopulib
make
```
But ideally, you should use this project as a template to create your own project. To do so,
please refer to the [CMakeTemplateRenamer](https://github.com/MattClarkson/CMakeTemplateRenamer)
which will show you how to clone this repository, and rename all the variables to names of your choice.
Then you would simply build your new project, using cmake, as shown above.


Further Build Instructions
--------------------------

This project can be configured, using CMake, to build against Eigen, Boost, OpenCV, glog, gflags, 
VTK, PCL and [ArrayFire](https://arrayfire.com/). These were chosen as examples of how to use CMake, and some common
C++ projects. These dependencies are optional, and this project will compile without them.

Furthermore, these dependencies can be downloaded and built,
or the user can specify directories of previously compiled
libraries.

To download and build dependencies, use CMake to ensure:

  * BUILD_SUPERBUILD:BOOL=ON

where ```ON``` is the default. Then to build any of Eigen, Boost or OpenCV etc., use CMake to set:

  * BUILD_Eigen:BOOL=ON|OFF
  * BUILD_Boost:BOOL=ON|OFF
  * BUILD_OpenCV:BOOL=ON|OFF

and so on. If you set BUILD_SUPERBUILD=OFF, and these BUILD_whatever variables are on, then CMake will just try finding
locally installed versions rather then downloading them.

To switch between static/dynamic linking, use CMake to set:

  * BUILD_SHARED_LIBS:BOOL=ON|OFF

To switch between Debug and Release mode, use CMake to set:

  * CMAKE_BUILD_TYPE:STRING=Debug|Release

Note: Only Debug and Release are supported.

Note: your host system is very likely to have a version of Boost that
is different to the one provided here. So if you want to use Boost,
you should probably try and use the one provided by this SuperBuild.


Caveats
-------

With all of those above build options, it is worth stressing:

 * You will still be required to write CMake code. This project is only to provide an EXAMPLE.
 * If you are building any of the dependencies, you would need to ensure the correct CMake flags are set to a reasonable default in order to compile those dependencies.
 * If you are testing on Travis or Appveyor, you need to configure your build to meet the required time limits, or else pay for more time.
 * If you are building Python Wheels, using the manylinux docker images, you may well need to provide your own docker image with dependencies pre-installed, or you setup correct CMake defaults such that the SuperBuild can compile the dependencies withing your docker image.

So, once more: This project is just to provide an example.
Its a template from which you can draw your own inspiration from.

 
Windows Users
-------------

If you build the project with shared libraries (BUILD_SHARED_LIBS:BOOL=ON)
then after the SuperBuild has successfully completed, you should look for the batch file
```StartVS_Debug.bat``` or ```StartVS_Release.bat``` in the ```PYVINECOPULIB-build``` folder.
This sets the path before launching Visual Studio, so that when you come to run your
application or unit tests within Visual Studio, the dynamically
loaded libraries are found at run time.


Running a Python Module
-----------------------

This project can be used to build Python extensions, using either Boost.Python
or PyBind11.

* Clone pyvinecopulib, using ```--recursive```.
* Use CMake to set BUILD_Python_Boost or BUILD_Python_PyBind to ON. These are mutually exclusive.
* Run a C++ build, as shown above.
* Set ```PYTHONPATH``` to the directory containing your C++ Python extension.

Code examples are in ```Code/PythonBoost``` or ```Code/PythonPybind```. 
So using the ```Code/PythonBoost``` as an example, once ```PYTHONPATH``` is set so
that you can pick up your compiled module, you would then be able to:
```
import  as mp
mp.my_first_add_function(1,6)
```

Deploying wheels is a difficult process. To help here, we have been 
inspired and assisted by [Matthew Brett's multibuild](https://github.com/matthew-brett/multibuild), 
so you would configure .travis.yml and appveyor.yml to generate your own Python Wheels.
This project itself is deployed to [PyPi](pypi.org) but obviously only includes the above
function for illustrative purposes.

Tested On
---------

 * Windows - Windows 8, VS2013, CMake 3.6.3, Qt 5.4.2
 * Linux - Centos 7, g++ 4.8.5, CMake 3.5.1, Qt 5.6.2
 * Mac - OSX 10.10.5, clang 6.0, CMake 3.9.4, Qt 5.6.2
 * Also refer to .travis.yml and appveyor.yml for other combinations

Minimum CMake version is 3.5. If you turn GUI options on, then the minimum Qt is version 5. 
Qt4 is not supported and not planned to be supported. If you are using VTK you should 
try a Qt version >= 5.5.0 to take advantage of the new OpenGL2 backend. 
If you need PyBind11, you need at least C++11, so on Windows you should have at least VS2015.


A Note on Packaging and Deployment
----------------------------------

There are many different issues and combinations to test when packaging an application.
For example:

 * System: Windows / Linux / Mac.
 * Linkage: Shared libraries / Static libraries.
 * Executable style: Command line applications / GUI applications / Command line applications bundled together with GUI applications.
 * With or without python modules.
 * With or without Unity modules.
 * The developer runs ```make install``` to install it in a specific directory, linked against known libraries, in known locations.
 * The developer runs ```make package``` to make a relocatable and distributable bundle so that another user can install it anywhere on their system.

It would take too long to document all of them and all of the issues involved. So, this project suggests
some simple starting points, and recommendations.

| Use Case | Important CMake settings | Workflow |
| -------- | ------------------------ | -------- |
| Command Line Apps | leave GUIs OFF, BUILD_SHARED_LIBS=OFF, CMAKE_INSTALL_PREFIX=/path/to/install/to | make, make install |
| C++ Gui Apps | Learn to build Qt5 first, BUILD_SHARED_LIBS=ON, build against your compiled Qt | make, make package | 
| Python Modules | must clone repo with --recursive, leave GUIs OFF, BUILD_SHARED_LIBS=OFF | make, test locally by setting PYTHONPATH, configure CI builds (see appveyor.yml and .travis.yml) to make wheels, upload to PyPi |
| Unity Modules |  leave GUIs OFF, BUILD_SHARED_LIBS=ON | make, then point Unity at the module. |

If you change your use-case, you must do a full clean build. That means completely delete the build folder, not just do a ```make clean```.


Preferred Branching Workflow for Contributions.
-----------------------------------------------

We welcome contributions to this project. Please use the following workflow.

 1. Raise issue in this project's Github Issue Tracker.
 2. Fork repository.
 3. Create a feature branch called ```<issue-number>-<some-short-description>```
    replacing ```<issue-number>``` with the Github issue number
    and ```<some-short-description>``` with your description of the thing you are implementing.
 4. Code on that branch.
 5. Push to your remote when ready.
 6. Create pull request.
 7. We will review code, suggest and required changes and merge to master when it is ready.


Examples Generated With This Template
-------------------------------------

 * [scikit-surgeryopencvcpp](https://github.com/UCL/scikit-surgeryopencvcpp) - image guided surgery functions using OpenCV, wrapped in Python.