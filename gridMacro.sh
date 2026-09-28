rm -rf 000
muse tarball > tarfile_$1.txt
tarname=$(cut -c 10- tarfile_$1.txt)
setup mu2etools
setup mu2efiletools
generate_fcl --dsconf=dummy       \
             --inputs=TrackerArt_09_21_DX_DT_Planev7_ZMinus.txt \
             --merge-factor=40        \
             --auto-description \
             --embed Tutorial/Alignment/fcl/sixpanelZ.fcl
tar czf $1.tar.bz2 000
#mkdir -p /pnfs/mu2e/scratch/users/$USER/fcl if never before
cp $1.tar.bz2 /pnfs/mu2e/scratch/users/$USER/fcl/
setup mu2egrid

mu2eprodsys --code=$tarname \
            --fcllist=/pnfs/mu2e/scratch/users/$USER/fcl/$1.tar.bz2  \
            --memory=10000MB \
            --expected-lifetime=24h \
            --xrootd \
            --dsconf=MDC2020v_perfect_v1_0
