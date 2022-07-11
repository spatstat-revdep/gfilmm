# gfilmm 2.0.5

Removed the 'rgr' dependency.


# gfilmm 2.0.4

Merged a pull request.


# gfilmm 2.0.3

The CRAN check shows an error. It is due to the update of the 'spatstat' package.


# gfilmm 2.0.2

The CRAN check on the Solaris platform still detect an error. This error does 
not occur on the Solaris platform of R-hub. I also checked my code with 
Valgrind, and no issue occured. So I am rather lost. Perhaps the issue is 
related to OpenMP (the package uses OpenMP for parallel computations). So, in 
this version, I force `nthreads = 1` if the platform is Solaris. I don't have 
any other idea, unfortunately.


# gfilmm 2.0.1

## Release summary

The CRAN checks of the previous version detected an error in the C++ code. I've 
found an error and I've fixed it, hopefully this error was the cause. I also 
fixed a NOTE ("all imports must be used").

## Test environments

* ubuntu 18.04, R 3.6.3
* win-builder (devel & release)
* r-hub

## R CMD check results

OK
