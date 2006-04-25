#!/bin/sh

subjdir=$1
slice=`expr $SGE_TASK_ID - 1`
slicezp=`${FSLDIR}/bin/zeropad $slice 4`

#$ -cwd -q long.q
#$ -S /bin/sh
#$ -V -N bedpostX
#$ -m a

${FSLDIR}/bin/xfibres\
 --data=$subjdir/data_slice_$slicezp\
 --mask=$subjdir/nodif_brain_mask_slice_$slicezp\
 -b $subjdir/bvals -r $subjdir/bvecs\
 --forcedir --logdir=$subjdir.bedpost/diff_slices/data_slice_$slicezp\
 --nj=1000 --bi=600 --se=20 --upe=24 --nfibres=2 > $subjdir.bedpost/logs/log$slicezp  && echo Done

