#!/bin/sh

#$ -cwd
#$ -q short.q
#$ -S /bin/sh
#$ -V
#$ -N bpx_postproc
#$ -m ae

subjdir=$1

numfib=`${FSLDIR}/bin/imglob -oneperimage ${subjdir}.bedpostX/diff_slices/data_slice_0000/f*samples | wc -w | awk '{print $1}'`

for ((fib=1; fib<=$numfib; fib++));do
    ${FSLDIR}/bin/avwmerge -z ${subjdir}.bedpostX/merged_th${fib}samples `${FSLDIR}/bin/imglob -oneperimage ${subjdir}.bedpostX/diff_slices/data_slice_*/th${fib}samples`
    ${FSLDIR}/bin/avwmerge -z ${subjdir}.bedpostX/merged_ph${fib}samples `${FSLDIR}/bin/imglob -oneperimage ${subjdir}.bedpostX/diff_slices/data_slice_*/ph${fib}samples`
    ${FSLDIR}/bin/avwmerge -z ${subjdir}.bedpostX/merged_f${fib}samples  `${FSLDIR}/bin/imglob -oneperimage ${subjdir}.bedpostX/diff_slices/data_slice_*/f${fib}samples`
    ${FSLDIR}/bin/avwmaths ${subjdir}.bedpostX/merged_th${fib}samples -Tmean ${subjdir}.bedpostX/mean_th${fib}samples
    ${FSLDIR}/bin/avwmaths ${subjdir}.bedpostX/merged_ph${fib}samples -Tmean ${subjdir}.bedpostX/mean_ph${fib}samples
    ${FSLDIR}/bin/avwmaths ${subjdir}.bedpostX/merged_f${fib}samples -Tmean ${subjdir}.bedpostX/mean_f${fib}samples

    ${FSLDIR}/bin/make_dyadic_vectors ${subjdir}.bedpostX/merged_th${fib}samples ${subjdir}.bedpostX/merged_ph${fib}samples ${subjdir}.bedpostX/nodif_brain_mask ${subjdir}.bedpostX/dyads${fib}
done

echo Removing intermediate files

if [ `imtest ${subjdir}.bedpostX/merged_th1samples` -eq 1 ];then
  if [ `imtest ${subjdir}.bedpostX/merged_ph1samples` -eq 1 ];then
    if [ `imtest ${subjdir}.bedpostX/merged_f1samples` -eq 1 ];then
      rm -rf ${subjdir}.bedpostX/diff_slices
      rm -f ${subjdir}/data_slice_*
      rm -f ${subjdir}/nodif_brain_mask_slice_*
    fi
  fi
fi

echo Creating identity xfm

xfmdir=${subjdir}.bedpostX/xfms
echo 1 0 0 0 > ${xfmdir}/eye.mat
echo 0 1 0 0 >> ${xfmdir}/eye.mat
echo 0 0 1 0 >> ${xfmdir}/eye.mat
echo 0 0 0 1 >> ${xfmdir}/eye.mat

echo Done
