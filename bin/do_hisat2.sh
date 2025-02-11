#!/bin/bash
i=1
hisat2_dir=$1
while true; do
    script1="$hisat2_dir/hisat2_ath${i}.sh"
    # 检查文件是否存在
    if [ -f "$script1" ]; then
        echo "sh $script1" 
    else
        break
    fi
    i=$((i + 1))
done > $hisat2_dir/all.hisat2.script.shell

cd $hisat2_dir
echo "/home/Project/ATGC_files/Pipeline/slurm/slurm_Duty.pl --interval 30 --maxjob 500 --convert no  --lines 1 --partition node001 --reslurm --mem 10G --cpu 2 --jobprefix work all.hisat2.script.shell" > qsub.sh
sh qsub.sh