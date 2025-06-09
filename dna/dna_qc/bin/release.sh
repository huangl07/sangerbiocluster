#!/bin/sh

mkdir -p $1/data_release/raw_data  $1/data_release/clean_data

mv  $1/01.CleanData/*clean*  $1/data_release/clean_data/

mv  $1/01.CleanData/*raw*  $1/data_release/raw_data/

#zip -r  $1/data_release/report.zip   $1/*report

mkdir -p $1/data_release/Report

mv  $1/WGS测序report/WGS测序report.html $1/data_release/Report/

mv $1/01.CleanData $1/data_release/Report/01.QC_Fig

rm -rf $1/data_release/Report/01.QC_Fig/*.list

rm -rf $1/data_release/Report/01.QC_Fig/qc/*.stat

mv $1/02.QCreport/qc-report.xls $1/data_release/Report/

rm -rf $1/02.QCreport/

rm -rf $1/WGS测序report/

cp $(dirname "$0")/readme/结果说明文档.txt $1/data_release/Report/
