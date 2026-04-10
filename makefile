# Single makefile for SlipOpt: builds nufft1d, nufft2d, rotgrid_r2014a
# 

OUTDIR=matlab

FC=/opt/homebrew/bin/flang
CC=/opt/homebrew/opt/llvm/bin/clang
CXX=/opt/homebrew/opt/llvm/bin/clang++

MATLAB_ROOT=/Applications/MATLAB_R2025b.app
MATLAB_INC=$(MATLAB_ROOT)/extern/include
MATLAB_LIB=$(MATLAB_ROOT)/bin/maca64

MW ?= ~/mwrap/mwrap
MWFLAGS=-c99complex -mex

FFLAGS=-g -O2 -fPIC -fopenmp -w

FLANG_RT_DIR ?= $(shell ls -d /opt/homebrew/opt/flang/lib/clang/*/lib/darwin 2>/dev/null | sort -V | tail -n 1)
ifeq ($(strip $(FLANG_RT_DIR)),)
$(error Could not locate flang runtime. Set FLANG_RT_DIR=/path/to/.../lib/darwin)
endif

MEX_LIBS=-L/opt/homebrew/opt/libomp/lib -lomp -L$(FLANG_RT_DIR) -lflang_rt.runtime -L$(MATLAB_LIB) -lmx -lmex -lmat -lm

NUFFT1D_O=$(OUTDIR)/nufft1df90.o $(OUTDIR)/nufft1dvf90.o $(OUTDIR)/dirft1d.o $(OUTDIR)/dfft.o $(OUTDIR)/next235.o
NUFFT2D_O=$(OUTDIR)/nufft2df90.o $(OUTDIR)/dirft2d.o $(OUTDIR)/dfft.o $(OUTDIR)/next235.o
ROTGRID_O=$(OUTDIR)/rotproj_cmpl.o $(OUTDIR)/rotviarecur3.o $(OUTDIR)/dfft.o $(OUTDIR)/fftnext.o $(OUTDIR)/yrecursion.o $(OUTDIR)/rotproj_real.o $(OUTDIR)/rotviarecur3_real.o $(OUTDIR)/rotproj_cmpl2.o $(OUTDIR)/rotproj_real2.o $(OUTDIR)/rotproj_cmpl3.o $(OUTDIR)/rotmat_proj.o $(OUTDIR)/rotviarecur3f.o $(OUTDIR)/rotviarecur3f_real.o $(OUTDIR)/sphtrans.o $(OUTDIR)/sphtrans_real.o $(OUTDIR)/legeexps.o $(OUTDIR)/chebexps.o $(OUTDIR)/prini.o $(OUTDIR)/rotgrid.o $(OUTDIR)/shiftphase.o $(OUTDIR)/nufft1df90.o $(OUTDIR)/nufft2df90.o $(OUTDIR)/nufft1dvf90.o $(OUTDIR)/next235.o $(OUTDIR)/nufft1df90_kb.o $(OUTDIR)/nufft2df90_kb.o $(OUTDIR)/i0eva.o $(OUTDIR)/rotgrid_aux.o $(OUTDIR)/gaussq.o
STOKES_O=$(OUTDIR)/sht_mod.o $(OUTDIR)/geometry_mod.o $(OUTDIR)/stokes_mod.o $(OUTDIR)/geometry_api_mex.o $(OUTDIR)/stokes_api_mex.o
STOKES_DEP_O=$(OUTDIR)/legeexps.o $(OUTDIR)/chebexps.o $(OUTDIR)/dfft.o $(OUTDIR)/prini.o $(OUTDIR)/rotproj_cmpl.o $(OUTDIR)/fftnext.o $(OUTDIR)/yrecursion.o $(OUTDIR)/sphtrans.o

.PHONY: all prep nufft1d nufft2d rotgrid stokes clean distclean

all: prep nufft1d nufft2d rotgrid stokes

prep:
	mkdir -p $(OUTDIR)
	rm -f $(OUTDIR)/*.o $(OUTDIR)/*.mod $(OUTDIR)/*.c $(OUTDIR)/*.mexmaca64

nufft1d: $(NUFFT1D_O)
	cp external/rotgrid/nufft1d/*.m $(OUTDIR)/
	cp external/rotgrid/nufft1d/nufft1d.mw $(OUTDIR)/
	cd $(OUTDIR) && $(MW) $(MWFLAGS) nufft1d -mb nufft1d.mw
	$(MW) $(MWFLAGS) nufft1d -c $(OUTDIR)/nufft1d.c external/rotgrid/nufft1d/nufft1d.mw
	$(CC) -c -fPIC -DMATLAB_MEX_FILE -DMWF77_UNDERSCORE1 -DMATLAB_DEFAULT_RELEASE=R2025a -I$(MATLAB_INC) $(OUTDIR)/nufft1d.c -o $(OUTDIR)/nufft1d.o
	$(CXX) -shared $(OUTDIR)/nufft1d.o $(NUFFT1D_O) $(MEX_LIBS) -o $(OUTDIR)/nufft1d.mexmaca64

nufft2d: $(NUFFT2D_O)
	cp external/rotgrid/nufft2d/*.m $(OUTDIR)/
	cp external/rotgrid/nufft2d/nufft2d.mw $(OUTDIR)/
	cd $(OUTDIR) && $(MW) $(MWFLAGS) nufft2d -mb nufft2d.mw
	$(MW) $(MWFLAGS) nufft2d -c $(OUTDIR)/nufft2d.c external/rotgrid/nufft2d/nufft2d.mw
	$(CC) -c -fPIC -DMATLAB_MEX_FILE -DMWF77_UNDERSCORE1 -DMATLAB_DEFAULT_RELEASE=R2025a -I$(MATLAB_INC) $(OUTDIR)/nufft2d.c -o $(OUTDIR)/nufft2d.o
	$(CXX) -shared $(OUTDIR)/nufft2d.o $(NUFFT2D_O) $(MEX_LIBS) -o $(OUTDIR)/nufft2d.mexmaca64

rotgrid: $(ROTGRID_O)
	$(MW) $(MWFLAGS) rotgrid_r2014a -c $(OUTDIR)/rotgrid_r2014a.c external/rotgrid/matlab/rotgrid_r2014a.mw
	$(CC) -c -fPIC -DMATLAB_MEX_FILE -DMWF77_UNDERSCORE1 -DMATLAB_DEFAULT_RELEASE=R2025a -I$(MATLAB_INC) $(OUTDIR)/rotgrid_r2014a.c -o $(OUTDIR)/rotgrid_r2014a.o
	$(CXX) -shared $(OUTDIR)/rotgrid_r2014a.o $(ROTGRID_O) $(MEX_LIBS) -o $(OUTDIR)/rotgrid_r2014a.mexmaca64

stokes: $(STOKES_O) $(STOKES_DEP_O)
	cd $(OUTDIR) && $(MW) $(MWFLAGS) stokes_mex -mb ../matlab/stokes.mw
	$(MW) $(MWFLAGS) stokes_mex -c matlab/stokes_mex.c matlab/stokes.mw
	$(CC) -c -fPIC -DMATLAB_MEX_FILE -DMWF77_UNDERSCORE1 -DMATLAB_DEFAULT_RELEASE=R2025a -I$(MATLAB_INC) matlab/stokes_mex.c -o $(OUTDIR)/stokes_mex.o
	$(CXX) -shared $(OUTDIR)/stokes_mex.o $(STOKES_O) $(STOKES_DEP_O) $(MEX_LIBS) -o matlab/stokes_mex.mexmaca64

$(OUTDIR)/sht_mod.o: src/sht_mod.f90
	$(FC) -c $(FFLAGS) -J$(OUTDIR) $< -o $@

$(OUTDIR)/stokes_mod.o: src/stokes_mod.f90 $(OUTDIR)/sht_mod.o
	$(FC) -c $(FFLAGS) -I$(OUTDIR) -J$(OUTDIR) $< -o $@

$(OUTDIR)/geometry_mod.o: src/geometry_mod.f90 $(OUTDIR)/sht_mod.o
	$(FC) -c $(FFLAGS) -I$(OUTDIR) -J$(OUTDIR) $< -o $@

$(OUTDIR)/geometry_api_mex.o: src/geometry_mex.f90 $(OUTDIR)/geometry_mod.o
	$(FC) -c $(FFLAGS) -I$(OUTDIR) $< -o $@

$(OUTDIR)/stokes_api_mex.o: src/stokes_mex.f90 $(OUTDIR)/stokes_mod.o
	$(FC) -c $(FFLAGS) -I$(OUTDIR) $< -o $@

$(OUTDIR)/%.o: external/rotgrid/src/%.f
	$(FC) -c $(FFLAGS) $< -o $@

$(OUTDIR)/%.o: external/rotgrid/nufft1d/%.f
	$(FC) -c $(FFLAGS) $< -o $@

$(OUTDIR)/%.o: external/rotgrid/nufft2d/%.f
	$(FC) -c $(FFLAGS) $< -o $@

clean:
	rm -f $(OUTDIR)/*.o $(OUTDIR)/*.mod $(OUTDIR)/*.mexmaca64 matlab/stokes_mex.mexmaca64

distclean:
	rm -rf $(OUTDIR)
