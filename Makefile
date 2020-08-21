include $(FSLCONFDIR)/default.mk

PROJNAME = fdt
SCRIPTS  = eddy_correct zeropad maskdyads probtrack fdt_rotate_bvecs \
           select_dwi_vols bedpost bedpostx bedpostx_postproc.sh \
           bedpostx_preproc.sh bedpostx_single_slice.sh \
           bedpostx_datacheck
FSCRIPTS = correct_and_average ocmr_preproc
XFILES   = dtifit ccops medianfilter make_dyadic_vectors vecreg xfibres \
           probtrackx pvmfit dtigen eddy_combine
FXFILES  = reord_OM sausages replacevols fdt_matrix_ops indexer \
           rearrange xfibres_pred
SOFILES  =
RUNTCLS  = Fdt
LIBS     = -lfsl-warpfns -lfsl-basisfield -lfsl-meshclass \
           -lfsl-newimage -lfsl-miscmaths -lfsl-NewNifti \
           -lfsl-utils -lfsl-znz -lfsl-cprob
CUDALIBS = -lcurand -lcudart -lcuda

ifeq ($(FDT_COMPILE_GPU), 1)
	XFILES  += merge_parts_gpu xfibres_gpu CUDA/split_parts_gpu
	SCRIPTS += CUDA/bedpostx_gpu CUDA/bedpostx_postproc_gpu.sh
    SOFILES += libfsl-bedpostx_cuda.so
endif


all: ${XFILES} ${FXFILES} ${SOFILES}

%: %.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

ccops: ccops.o ccopsOptions.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

probtrackx: probtrackx.o probtrackxOptions.o streamlines.o ptx_simple.o ptx_seedmask.o ptx_twomasks.o ptx_nmasks.o ptx_meshmask.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

dtifit: dtifit.o dtifitOptions.o diffmodels.o Bingham_Watson_approx.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

xfibres: xfibres.o xfibresoptions.o diffmodels.o Bingham_Watson_approx.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

xfibres_2: xfibres_2.o xfibresoptions.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

kurtosis: kurtosis.o dtifitOptions.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

pvmfit: pvmfit.o pvmfitOptions.o diffmodels.o Bingham_Watson_approx.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

basgen: basgen.o diffmodels.o Bingham_Watson_approx.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

rubix: rubix.o diffmodels.o rubixvox.o rubixoptions.o Bingham_Watson_approx.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

libfsl-bedpostx_cuda.so:
	${NVCC} --shared \
	  ${CUDACXXFLAGS} -ICUDA -ICUDA/options \
      -o $@ \
      CUDA/init_gpu.cu CUDA/samples.cu CUDA/diffmodels.cu \
      CUDA/runmcmc.cu CUDA/xfibres_gpu.cu ${CUDALDFLAGS}

merge_parts_gpu: merge_parts_gpu.o xfibresoptions.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

CUDA/split_parts_gpu: CUDA/split_parts_gpu.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

XFIBRES_OBJS = xfibres_gpu.o xfibresoptions.o diffmodels.o \
               Bingham_Watson_approx.o

xfibres_gpu: libfsl-bedpostx_cuda.so ${XFIBRES_OBJS}
	${CXX} ${CXXFLAGS} -o $@ ${XFIBRES_OBJS} \
        ${LDFLAGS} -lfsl-bedpostx_cuda ${CUDALDFLAGS}
