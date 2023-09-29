THISDIR=$PWD
MYBUILD=../../../build

cd $MYBUILD
#make -j10 install

cd $THISDIR
echo $THISDIR



cd g4ms_output
rm shortlog.log
echo 'ciao' > shortlog.log

#awk '/Gun                               INFO/,EOF { print $0 }' g4ms_*.log > shortlog.log


awk '/Gun                               INFO/,EOF { print $0 }' g4ms_A_N.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_A_E.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_A_S.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_A_W.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_B_N.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_B_E.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_B_S.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_B_W.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_C_N.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_C_E.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_C_S.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_C_W.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_D_N.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_D_E.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_D_S.log >> shortlog.log
awk '/Gun                               INFO/,EOF { print $0 }' g4ms_D_W.log >> shortlog.log


awk '/Gun/ || /Air/ || /AlBeMet/ || /Au/ || /Cu/ || /LiquidNDecane/ || /Radiation/ || /Length/ || /Layer/ {print $0}' shortlog.log  
exit

checkOverlaps -c standalone_CADbp.xml

#sample points (x,y,z in mm):
#
# A: 0,0,0       IP
# B: 0,0,125     central cooling manifold
# C: 0,0,310     long  trap cooling manifold (+10mm)
# D: 0,0,630     short trap cooling manifold (+10mm)

#sample directions: N-E-S-W (compass)
mkdir g4ms_output

g4MaterialScan -c standalone_CADbp.xml -d 0,1,0  -p 0,0,0 > g4ms_output/g4ms_A_N.log
g4MaterialScan -c standalone_CADbp.xml -d 1,0,0  -p 0,0,0 > g4ms_output/g4ms_A_E.log
g4MaterialScan -c standalone_CADbp.xml -d 0,-1,0 -p 0,0,0 > g4ms_output/g4ms_A_S.log
g4MaterialScan -c standalone_CADbp.xml -d -1,0,0 -p 0,0,0 > g4ms_output/g4ms_A_W.log

g4MaterialScan -c standalone_CADbp.xml -d 0,1,0  -p 0,0,125 > g4ms_output/g4ms_B_N.log
g4MaterialScan -c standalone_CADbp.xml -d 1,0,0  -p 0,0,125 > g4ms_output/g4ms_B_E.log
g4MaterialScan -c standalone_CADbp.xml -d 0,-1,0 -p 0,0,125 > g4ms_output/g4ms_B_S.log
g4MaterialScan -c standalone_CADbp.xml -d -1,0,0 -p 0,0,125 > g4ms_output/g4ms_B_W.log

g4MaterialScan -c standalone_CADbp.xml -d 0,1,0  -p 0,0,310 > g4ms_output/g4ms_C_N.log
g4MaterialScan -c standalone_CADbp.xml -d 1,0,0  -p 0,0,310 > g4ms_output/g4ms_C_E.log
g4MaterialScan -c standalone_CADbp.xml -d 0,-1,0 -p 0,0,310 > g4ms_output/g4ms_C_S.log
g4MaterialScan -c standalone_CADbp.xml -d -1,0,0 -p 0,0,310 > g4ms_output/g4ms_C_W.log

g4MaterialScan -c standalone_CADbp.xml -d 0,1,0  -p 0,0,630 > g4ms_output/g4ms_D_N.log
g4MaterialScan -c standalone_CADbp.xml -d 1,0,0  -p 0,0,630 > g4ms_output/g4ms_D_E.log
g4MaterialScan -c standalone_CADbp.xml -d 0,-1,0 -p 0,0,630 > g4ms_output/g4ms_D_S.log
g4MaterialScan -c standalone_CADbp.xml -d -1,0,0 -p 0,0,630 > g4ms_output/g4ms_D_W.log


