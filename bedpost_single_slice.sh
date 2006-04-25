#!/bin/sh

subjdir=$1
slice=`expr $SGE_TASK_ID - 1`
slicezp=`${FSLDIR}/bin/zeropad $slice 4`

#$ -cwd -q long.q
#$ -S /bin/sh
#$ -V -N bedpost
#$ -m a

${FSLDIR}/bin/diff_pvm\
 --data=$subjdir/data_slice_$slicezp\
 --mask=$subjdir/nodif_brain_mask_slice_$slicezp\
 -b $subjdir/bvals -r $subjdir/bvecs\
 --forcedir --logdir=$subjdir.bedpost/diff_slices/data_slice_$slicezp\
 --nj=1300 --bi=300 --se=20 > $subjdir.bedpost/logs/log$slicezp && echo Done
