clear
for f in *000.launch; do 
    MOL=$(echo $f | cut -d"_" -f1)
    NUM=$(ls ${MOL}_pbe_?????_*.out | wc -l)
    if [ ! "$NUM" == "0" ]; then 
        for OUTPUT in ${MOL}_pbe_?????_*.out; do
            JOB=$(echo $OUTPUT | cut -d"." -f1)
            if [ "$(normaltermination_mrchem.py $OUTPUT)" == "False" ]; then 
                ORBDIR=/global/work/ambr/benchmark_orbitals/0001_highprec/${MOL}_pbe_000
                LAUNCHFILE=${JOB}.launch
                
                echo $LAUNCHFILE
                sed -i s/"--cpus-per-task=16"/"--cpus-per-task=8"/ $LAUNCHFILE
                sbatch $LAUNCHFILE
            fi
        done
    fi
done 2> /dev/null
