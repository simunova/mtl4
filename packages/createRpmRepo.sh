mkdir -p mtl_rpm_repo/main
cp standard/*.rpm mtl_rpm_repo/main/
cp examples/*.rpm mtl_rpm_repo/main/
cp tests/*.rpm mtl_rpm_repo/main/
cp all/*.rpm mtl_rpm_repo/main/
cd mtl_rpm_repo/main
createrepo .
gpg -a --detach-sign repodata/repomd.xml
gpg -a --export > repodata/repomd.xml.key
cd ../../
