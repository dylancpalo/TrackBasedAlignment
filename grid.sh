rm -rf 000
muse tarball PassN/TrackerCalib_Palo/TrkPanelMap.txt PassN/TrackerCalib_Palo/OtherTiming.txt PassN/TrackerCalib_Palo/PashaPanelTiming_v4.txt PassN/TrackerCalib_Palo/dukeOnly.txt PassN/TrackerCalib_Palo/PanelAlignment_v3.txt PassN/TrackerCalib_Palo/PlaneAlignv7_ZMinus.txt > tarfile_$1.txt
tarname=$(cut -c 10- tarfile_$1.txt)
setup mu2etools
setup mu2efiletools
generate_fcl --dsconf=dummy       \
             --inputs=Tracker124.txt \
             --merge-factor=1        \
             --auto-description \
             --embed PassN/TrackerCalib_Palo/TrkRecoFromFragments.fcl
tar czf $1.tar.bz2 000
#mkdir -p /pnfs/mu2e/scratch/users/$USER/fcl if never before
cp $1.tar.bz2 /pnfs/mu2e/scratch/users/$USER/fcl/
setup mu2egrid

mu2eprodsys --code=$tarname \
            --fcllist=/pnfs/mu2e/scratch/users/$USER/fcl/$1.tar.bz2  \
            --memory=10000MB \
            --expected-lifetime=12h \
            --xrootd \
            --dsconf=MDC2020v_perfect_v1_0
