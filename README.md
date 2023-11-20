# MPI-Network-Resource-optimization


To run the code:
We used the A DMC server for our approach. However, it should run on all other servers (picocluster and Jetson).

Since we used MPI, the first step is to run:
module load openmpi

After creating the MPI environment, you can create a gnu file using this command:
mpic++ project.cpp -o project_gnu

This will create a project_gnu file. Now we need to create a shell file to execute the program. We created a shell file named project_gnu.sh and wrote:
#!/bin/bash
module load openmpi
srun --mpi=pmi2 ./project_gnu 100 13 1000.0 0.1 >project.txt


After that, we gave permission to the .sh file using
chmod 755 project_gnu.sh

Then we ran our code in DMC server using
run_script_mpi project_gnu.sh

This will ask for some input values including

Enter Queue Name (default <cr>: small) class
Enter number of processor cores (default <cr>: 1 ) 1/2/3/4/8/16
Enter Time Limit (default <cr>: 12:00:00 HH:MM:SS) 
Enter memory limit (default <cr>: 16gb ) 
Enter a name for your job (default: projectgnushSCRIPT)
Run on; broadwell, skylake, milan (default: dmc)
