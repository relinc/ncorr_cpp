#NCORR for Windows
This build is for Windows. 
To build, create a build directory, run the following command in the power shell from that directory (once cmake and visual studio are installed).
```
cmake -G "Visual Studio 15 2017 Win64" ..
```

Then open the .sln in Visual Studio, right click on the build all project select "Run Build".
Then build the install option as well.

Default settings:
Lagrangian strain is chosen for the output strain as it can be easily transformed to true strain. Lagrangian strain is defined with respect to the initial length of the the object, which is the usual for mechanics.

Eulerian strain is chosen for video as it allows for easy display of a deformed strain field. This allows for easier visualization of how the material is deforming.

Scale factor is the factor that displacement results are downsized. A larger factor gives faster, less accurate results.

Interpolation, cubic keys interpolation is defaulted to as it saves time on computation without a large degradation in quality. Biquintic interpolation is more accurate, but takes longer to compute.

Threads is the number of threads the correlation job is broken into. The more threads, the faster the job completes, however, less resources are left for using the computer while computation is being done.

Subregion is defaulted to circle. In the future, rectangular subregions may be offered. The subregion is the shape small area tracked by the correlation.

Radius sub is the radius of the subregion being tracked. This must correspond to have a few identifiable markers in the radius.



