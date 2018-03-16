#!/bin/bash

# Note: does not work with Ensembl release < 82

RELEASE=$1
URL_FTP=ftp://ftp.ensembl.org/pub/release-${RELEASE}
FILE_GTF=Homo_sapiens.GRCh38.${RELEASE}.chr.gtf.gz

wget ${URL_FTP}/gtf/homo_sapiens/${FILE_GTF}
echo Rscript scripts/prepare_anno_ensembl.R ${FILE_GTF} ${RELEASE}

rm -rf ${FILE_GTF}
