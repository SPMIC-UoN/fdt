#!/bin/sh

#$ -cwd 
#$ -q long.q
#$ -S /bin/sh
#$ -V
#$ -N bedpostx
#$ -m a
#$ -o $1.bedpostX/logs -e $1.bedpostX/logs

subjdir=$1
nfibres=$2
fudge=$3
bi=$4
slice=$5

slicezp=`${FSLDIR}/bin/zeropad $slice 4`

${FSLDIR}/bin/xfibres\
 --data=$subjdir/data_slice_$slicezp\
 --mask=$subjdir/nodif_brain_mask_slice_$slicezp\
 -b $subjdir/bvals -r $subjdir/bvecs\
 --forcedir --logdir=$subjdir.bedpostX/diff_slices/data_slice_$slicezp\
 --fudge=$fudge --nj=5000 --bi=$bi --se=100 --upe=24 --nfibres=$nfibres > $subjdir.bedpostX/logs/log$slicezp  && echo Done

