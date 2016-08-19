
for mach in $(seq 0.8 0.1 3.0)
do 
    ./mk-config-dimless.tcl flash.par.stability.ttcl $mach \
        > /mnt/sshfs/markert-srv/FLASH/stirturb/es/2016-07-11/stability_region/flash.par/$mach
done
