#!/bin/sh
rm -rf downloadtable
mkdir -p downloadtable

SUFFICES=".tar.gz;TGZ .tar.bz2;TBZ2 .rpm;RPM .zip;ZIP .deb;DEB";
cp downloadtable_base.html downloadtable/downloadtable.html
for suffs in ${SUFFICES} ; do
	filetype=`echo ${suffs} | cut -d';' -f1`
	namesuff=`echo ${suffs} | cut -d';' -f2`
	echo "filetype: ${filetype}"
	echo "namesuff: ${namesuff}"

	cp standard/MTL-*${filetype} downloadtable/
	FILE_STANDARD=`ls standard/MTL-*${filetype} | sed 's:standard/::'`
	cp examples/MTL-examples*${filetype} downloadtable/
	FILE_EXAMPLES=`ls examples/MTL-examples-*${filetype} | sed 's:examples/::'`
	cp tests/MTL-tests*${filetype} downloadtable/
	FILE_TESTS=`ls tests/MTL-tests-*${filetype} | sed 's:tests/::'`
	cp all/MTL-all*${filetype} downloadtable/
	FILE_ALL=`ls all/MTL-all-*${filetype} | sed 's:all/::'`

	MD5_STANDARD=`md5sum standard/MTL-*${filetype} | sed 's/ .*$//'` ;
	MD5_EXAMPLES=`md5sum examples/MTL-examples-*${filetype} | sed 's/ .*$//'` ;
	MD5_TESTS=`md5sum tests/MTL-tests-*${filetype} | sed 's/ .*$//'`;
	MD5_ALL=`md5sum all/MTL-all*${filetype} | sed 's/ .*$//'`;

	SIZE_STANDARD=`ls -sh downloadtable/${FILE_STANDARD} | cut -d' ' -f1` ;
	SIZE_EXAMPLES=`ls -sh downloadtable/${FILE_EXAMPLES} | cut -d' ' -f1` ;
	SIZE_TESTS=`ls -sh downloadtable/${FILE_TESTS} | cut -d' ' -f1` ;
	SIZE_ALL=`ls -sh downloadtable/${FILE_ALL} | cut -d' ' -f1` ;

	sed -i "s:MTL4_MD5_${namesuff}:${MD5_STANDARD}:g" downloadtable/downloadtable.html
	sed -i "s:EXAMPLES_MD5_${namesuff}:${MD5_EXAMPLES}:g" downloadtable/downloadtable.html
	sed -i "s:TESTS_MD5_${namesuff}:${MD5_TESTS}:g" downloadtable/downloadtable.html
	sed -i "s:ALL_MD5_${namesuff}:${MD5_ALL}:g" downloadtable/downloadtable.html
	sed -i "s:MTL4_FILE_${namesuff}:${FILE_STANDARD}:g" downloadtable/downloadtable.html
	sed -i "s:EXAMPLES_FILE_${namesuff}:${FILE_EXAMPLES}:g" downloadtable/downloadtable.html
	sed -i "s:TESTS_FILE_${namesuff}:${FILE_TESTS}:g" downloadtable/downloadtable.html
	sed -i "s:ALL_FILE_${namesuff}:${FILE_ALL}:g" downloadtable/downloadtable.html
	sed -i "s:MTL4_SIZE_${namesuff}:${SIZE_STANDARD}:g" downloadtable/downloadtable.html
	sed -i "s:EXAMPLES_SIZE_${namesuff}:${SIZE_EXAMPLES}:g" downloadtable/downloadtable.html
	sed -i "s:TESTS_SIZE_${namesuff}:${SIZE_TESTS}:g" downloadtable/downloadtable.html
	sed -i "s:ALL_SIZE_${namesuff}:${SIZE_ALL}:g" downloadtable/downloadtable.html
done
exit 0
