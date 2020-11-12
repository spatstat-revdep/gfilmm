# gfilmm 1.0.1

* Removed the quadruple precision feature because of issues with the CRAN checks.
* The main function `gfilmm` now has an option `long`, to run the algorithm in 
long double precision instead of double precision.


# gfilmm 1.0.0

* Fixed two mistakes in the C++ code.
* Fixed a mistake in the R code; the results were wrong when the model had random effects with interaction.
* The main function `gfilmm` now has an option `precision`, to run the algorithm in long double precision or quadruple precision instead of double precision.
* It also has a new option seed, to set the seed of the C++ random engine.
* Updated the vignette.


# gfilmm 0.1.0

First release.
