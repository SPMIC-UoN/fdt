#!/bin/sh

subjdir=$1

#$ -S /bin/sh
#$ -V
#$ -N bp_preproc
#$ -m a

echo Copying files to bedpost directory
cp ${subjdir}/bvecs ${subjdir}/bvals ${subjdir}.bedpost
${FSLDIR}/bin/imcp ${subjdir}/nodif ${subjdir}/nodif_brain_mask ${subjdir}.bedpost
${FSLDIR}/bin/avwmaths\
 ${subjdir}/nodif\
 -mas ${subjdir}/nodif_brain_mask\
 ${subjdir}.bedpost/nodif_brain

${FSLDIR}/bin/avwslice ${subjdir}/data
${FSLDIR}/bin/avwslice ${subjdir}/nodif_brain_mask
echo Done