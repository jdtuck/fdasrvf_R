## Test environments
* local OS X install, R 3.2.2
* ubuntu 12.04 (on travis-ci), R 3.2.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs and 1 WARNINGs:

* Found the following file(s) containing GNU extensions:
  src/Makevars Portable Makefiles do not use GNU extensions such as +=, :=, $(shell), $(wildcard), ifeq ... endif. See section ‘Writing portable packages’ in the ‘Writing R Extensions’ manual.

  The dependent c code has a lot of source files and I am using a subdirectory
  structure for better organization

There was 1 NOTE:

* New maintainer:
  J. Derek Tucker <jdtuck@sandia.gov>
  Old maintainer(s):
  J. Derek Tucker <dtucker@stat.fsu.edu>

  I have a change in e-mail addresses

## Downstream dependencies
There are no downstream dependencies
