#!/bin/bash


dataarray=('mnist8m' 'brats17')
bwarray=(0.42 0.35) 
regarray=(1.0E-6 9.4E-8)
rankarray=(8196 16384)
rankarray=(8192)

#dataarray=('sphere4')
#bwarray=(0.56)
#regarray=(4.9E-6)
#rankarray=(512 1024 2048)

for ((r=0;r<${#rankarray[@]};++r)); do
    rank=${rankarray[r]}
    for ((i=0;i<${#dataarray[@]};++i)); do
    	cur_data=${dataarray[i]}
    	cur_bw=${bwarray[i]}
    	cur_reg=${regarray[i]}
    	echo "Setting up run with ${cur_data} at bw ${cur_bw} and reg ${cur_reg} for rank ${rank}"
    	
    	# Make sbatch runfile
    	echo "#!/bin/bash
#SBATCH -J rcd_${cur_data}_${rank}            # job name
#SBATCH -o rcd_${cur_data}_${rank}.o        # output and error file name (%j expands to jobID)
#SBATCH -N 1                # number of nodes requested
#SBATCH -n 1               # total number of mpi tasks requested
#SBATCH -p largemem512GB   
#SBATCH -t 4:30:00         # run time (hh:mm:ss) - 1.5 hours
# Slurm email notifications are now working on Lonestar 5 
#SBATCH --mail-user=sameer@ices.utexas.edu
#SBATCH --mail-type=end     # email me when the job finishes
#SBATCH --mail-type=begin     # email me when the job finishes


pwd
module load matlab/2019a
cd /work/03158/tharakan/research/rklr/matlab/ 


matlab -r \"set_local_env; run_rc_experiments('${cur_data}',${cur_bw},${cur_reg},${rank} ); quit\"


" > ${cur_data}.${rank}job
	
	
	sbatch ${cur_data}.${rank}job
    done
	
done

