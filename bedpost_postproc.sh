#!/bin/sh

subjdir=$1

#$ -S /bin/sh
#$ -V
#$ -N bp_postproc
#$ -m ae

echo Merging outputs into 4D files

${FSLDIR}/bin/avwmerge -z ${subjdir}.bedpost/merged_thsamples `${FSLDIR}/bin/imglob -oneperimage ${subjdir}.bedpost/diff_slices/data_slice_*/th_samples`
${FSLDIR}/bin/avwmerge -z ${subjdir}.bedpost/merged_phsamples `${FSLDIR}/bin/imglob -oneperimage ${subjdir}.bedpost/diff_slices/data_slice_*/ph_samples`
${FSLDIR}/bin/avwmerge -z ${subjdir}.bedpost/merged_fsamples  `${FSLDIR}/bin/imglob -oneperimage ${subjdir}.bedpost/diff_slices/data_slice_*/f_samples`
${FSLDIR}/bin/avwmaths ${subjdir}.bedpost/merged_thsamples -Tmean ${subjdir}.bedpost/mean_thsamples
${FSLDIR}/bin/avwmaths ${subjdir}.bedpost/merged_phsamples -Tmean ${subjdir}.bedpost/mean_phsamples
${FSLDIR}/bin/avwmaths ${subjdir}.bedpost/merged_fsamples -Tmean ${subjdir}.bedpost/mean_fsamples

${FSLDIR}/bin/make_dyadic_vectors ${subjdir}.bedpost/merged_thsamples ${subjdir}.bedpost/merged_phsamples ${subjdir}.bedpost/nodif_brain_mask ${subjdir}.bedpost/dyadic_vectors

echo Removing intermediate files

if [ `imtest ${subjdir}.bedpost/merged_thsamples` -eq 1 ];then
  if [ `imtest ${subjdir}.bedpost/merged_phsamples` -eq 1 ];then
    if [ `imtest ${subjdir}.bedpost/merged_fsamples` -eq 1 ];then
      rm -rf ${subjdir}.bedpost/diff_slices
      rm -f ${subjdir}/data_slice_*
      rm -f ${subjdir}/nodif_brain_mask_slice_*
    fi
  fi
fi

echo Creating identity xfm

xfmdir=${subjdir}.bedpost/xfms
echo 1 0 0 0 > ${xfmdir}/eye.mat
echo 0 1 0 0 >> ${xfmdir}/eye.mat
echo 0 0 1 0 >> ${xfmdir}/eye.mat
echo 0 0 0 1 >> ${xfmdir}/eye.mat

echo Done