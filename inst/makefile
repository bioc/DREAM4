default: quickbuild install

help:
	egrep "^#" makefile | sed "s/^#//"

# --- quickbuild: no vignette
#
quickbuild:
	(cd ../..; R CMD build --no-vignettes DREAM4)


# --- build
#
build:
	(cd ../..; R CMD build DREAM4)

# --- install
#
install: quickbuild
	(cd ../..; R CMD install DREAM4)

# --- check
#
check: build install
	(cd ../..; R CMD check --no-manual --no-vignettes --no-codoc --no-examples --no-manual  DREAM4)

# --- checkfull
#
checkfull: build install
	(cd ../..; R CMD check DREAM4)


# --- vanillaTest
# run all the unit tests, in a clean context
#

vanillaTest:  build install
	- rm vanillaTest.out
	R --vanilla < vanillaTest.R > vanillaTest.out 2>&1

# --- vt
# run all the unit tests, in a clean context
#

vt: vanillaTest



# --- checkvig
# check just the vignette
#

checkvig:
	(cd ../..; R CMD check --no-manual --no-codoc --no-tests --no-examples DREAM4)


# --- tangle
# extract the R code from the vignette file
#

tangle:
	(cd doc; R CMD Stangle $(PACKAGE_NAME).Rnw)



# --- sweave
# creates DREAM4.tex, runs all embedded examples
# run this before the pdf target
#
sweave: 
	(cd doc; R CMD Sweave $(PACKAGE_NAME).Rnw)

# --- pdf
# make and open DREAM4.pdf, the vignette 
#

pdf:  sweave
	(cd doc; rm *.pdf; /usr/texbin/latex --output-format=pdf  DREAM4.tex; open DREAM4.pdf)



