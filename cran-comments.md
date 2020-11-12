## Release summary

With the previous version, the CRAN checks returned an error on systems using 
the 'clang' compiler. This is because 'clang' has no support for the 'float128' 
numbers of the C++ Boost library. So I've removed the possibility to use these 
numbers.

## Test environments

* ubuntu 18.04, R 3.6.3
* win-builder (devel & release)
* r-hub

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
