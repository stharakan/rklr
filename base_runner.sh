#!/bin/bash

dataarray=('sphere4' 'mnist8m' 'brats17')

dataarray=( 'sphere4' )

for ((i=0;i<${#dataarray[@]};++i)); do
	cur_data=${dataarray[i]}
	echo "Setting up run with ${cur_data} "
	
	# Make sbatch runfile
	echo "#!/bin/bash
#SBATCH -J bw_${cur_data}_512            # job name
#SBATCH -o bw_${cur_data}_512.o        # output and error file name (%j expands to jobID)
#SBATCH -N 1                # number of nodes requested
#SBATCH -n 1               # total number of mpi tasks requested
#SBATCH -p normal 
##SBATCH -p largemem512GB
#SBATCH -t 1:40:00         # run time (hh:mm:ss) - 1.5 hours
# Slurm email notifications are now working on Lonestar 5 
#SBATCH --mail-user=sameer@ices.utexas.edu
#SBATCH --mail-type=end     # email me when the job finishes


pwd
module load matlab/2019a
cd /work/03158/tharakan/research/rklr/matlab/ 


matlab -r \"set_local_env; run_sigma_selection('${cur_data}',512 ); quit\"  

" > ${cur_data}.job
	
	
	sbatch ${cur_data}.job
	
done

