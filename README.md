Prerequsites:

GLib

> Specifically: `libglib2.0-dev`

Files:

> Narative story of process of measuring and modifying code
(start with this one) 
```
Analysis.pdf
```
> Source for modified, optimized code

```
jstate.h
main.c
thpool.c
thpool.h
```

> Folder of saved profiler, and Calgrind logs 

```
PerfTests
```
> Dump of output over coarse grid.  Used to validate accuracy of faster versions

```
smalltruth.txt

```

To build single threaded version:

> Clone repo.

> Default make target is non-printing "release" version taking step size on command line. Default is coarse undersampling 0.2.  Original value is 0.01.
	
	
> Build and test as described in document `Analysis.pdf`

For multi threaded version:

> Checkout branch called threadpool, build and test as described in Analysis document mentioned above. 
	
	


