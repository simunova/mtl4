#!/bin/sh

cp ../standard/MTL-*.tar.gz ../examples/MTL-examples*.tar.gz ../tests/MTL-tests*.tar.gz . ;
LONGNAME=`ls MTL-all-*.tar.gz | sed 's/\.gz//'` ;
ALL_DIRNAME=`echo ${LONGNAME} | sed 's/\.tar//'` ;

rm MTL-all*.tar.bz2 ;
rm MTL-all*.tar.gz ;

mkdir ${ALL_DIRNAME} ;
for d in $(ls MTL-*.tar.gz); do
tar -xf ${d} ;
rm -rf ${d} ;
curdirname=`echo ${d} | sed 's/\.tar\.gz//'` ;
cp -r ${curdirname}/* ${ALL_DIRNAME}/ ;
rm -rf ${curdirname} ;
done

tar -czf ${ALL_DIRNAME}.tar.gz ${ALL_DIRNAME}
tar -cjf ${ALL_DIRNAME}.tar.bz2 ${ALL_DIRNAME}
zip -q -r ${ALL_DIRNAME}.zip ${ALL_DIRNAME}
