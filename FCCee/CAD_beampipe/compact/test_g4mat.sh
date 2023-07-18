THISDIR=$PWD
MYBUILD=../../../build

cd $MYBUILD
make -j10 install

cd $THISDIR
echo $THISDIR

g4MaterialScan -c standalone_CADbp.xml -d 0.05,0,1.0
