[![CBox with smoke](https://github.com/moonlovelj/nori/blob/master/scenes/volume/cbox_smoke.png)](https://github.com/moonlovelj/nori/blob/master/scenes/volume/cbox_smoke.png)
[![vech](https://github.com/moonlovelj/nori/blob/master/scenes/pa5/veach_mi/veach_mis.png)](https://github.com/moonlovelj/nori/blob/master/scenes/pa5/veach_mi/veach_mis.png)
[![table](https://github.com/moonlovelj/nori/blob/master/scenes/pa5/table/table_mis.png)](https://github.com/moonlovelj/nori/blob/master/scenes/pa5/table/table_mis.png)
[![table rouch dielectic](https://github.com/moonlovelj/nori/blob/master/scenes/pa5/table/table_ajax_mis.png)](https://github.com/moonlovelj/nori/blob/master/scenes/pa5/table/table_ajax_mis.png)

## Nori Version 2
![Build status](https://github.com/wjakob/nori/workflows/Build/badge.svg)

Nori is a simple ray tracer written in C++. It runs on Windows, Linux, and
Mac OS and provides basic functionality that is required to complete the
assignments in the course Advanced Computer Graphics taught at EPFL.

### Course information and framework documentation

For access to course information including slides and reading material, visit the main [Advanced Computer Graphics website](https://rgl.epfl.ch/courses/ACG17). The Nori 2 framework and coding assignments will be described on the [Nori website](https://wjakob.github.io/nori).

### Note to researchers and students from other institutions

Last year's version of Nori including a full set of assignment descriptions can
be found in the following [archive](https://github.com/wjakob/nori-old).


### Known Issues
There is a known issue with the NanoGUI version that Nori uses: on Linux systems with an integrated Intel GPU, a bug in the Mesa graphics drivers causes the GUI to freeze on startup. A workaround is to temporarily switch to an older Mesa driver to run Nori. This can be done by running
```
export MESA_LOADER_DRIVER_OVERRIDE=i965
```
