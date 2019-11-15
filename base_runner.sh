

dataarray=('sphere4' 'brats14' 'mnist8m')
bwarray=(0.56 0.42 0.35) 


for ((i=0;i<${#dataarray[@]};++i)); do
	cur_data=${dataarray[i]}
	cur_bw=${bwarray[i]}
	echo "Setting up run with ${cur_data} at bw ${bw}"
	
	# Make sbatch runfile
	echo "#!/bin/bash
	#SBATCH -J reg_${cur_data}_1024            # job name
	#SBATCH -o reg_${cur_data}_1024.o        # output and error file name (%j expands to jobID)
	#SBATCH -N 1                # number of nodes requested
	#SBATCH -n 1               # total number of mpi tasks requested
	#SBATCH -p largemem512GB
	#SBATCH -t 0:40:00         # run time (hh:mm:ss) - 1.5 hours
	# Slurm email notifications are now working on Lonestar 5 
	#SBATCH --mail-user=sameer@ices.utexas.edu
	#SBATCH --mail-type=end     # email me when the job finishes
	
	
	pwd
	module load matlab/2019a
	cd /work/03158/tharakan/research/rklr/matlab/ 
	
	
	matlab -r \"set_local_env; run_lambda_selection('${cur_data}',1024,${cur_bw} ); quit\"  
	
	" > ${cur_data}.job
	
	
	sbatch ${cur_data}.job
	
done

