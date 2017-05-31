#!/usr/bin/env make -f
# Makefile called by {SPM}/external/Makefile
#
# 
# $Id: Makefile 1026 2013-02-05 19:15:21Z guillaume $

include ../../src/Makefile.var

ACIDMEX = \
	cfiles/spm_hist2_z_exp_polyval2.$(SUF) cfiles/DTI2EVEW.$(SUF) \
	cfiles/DTtoFA_2d.$(SUF) cfiles/DTtoFA_nachHBM2010.$(SUF) \
	Preprocessing/HySCO/n2ccScalarMexC.$(SUF) \
	Preprocessing/HySCO/FAIRkernel/splineInterMexC.$(SUF) \
	Preprocessing/HySCO/FAIRkernel/linearInterMexC.$(SUF) \
	Preprocessing/HySCO/getPartialBMexC.$(SUF)

all: $(ACIDMEX)

clean:

distclean: clean
	$(DEL) $(ACIDMEX)

install:

tarball: all
	$(TAR) -cf acid_mex.tar $(ACIDMEX)

cfiles/%.$(SUF) : cfiles/%.c
	$(MEX) $< -outdir $(dir $<) $(MEXEND)

Preprocessing/HySCO/%.$(SUF) : Preprocessing/HySCO/%.cpp
	$(MEX) $< -outdir $(dir $<) $(MEXEND)
