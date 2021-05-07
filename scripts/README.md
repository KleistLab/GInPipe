# How to run the test simulations pipeline

The testing pipeline can be executed using **call_evol_run.sh**. This is bash script. To run it, set the working directory and results directory variables in the top of the file:

```
WORK_DIR="path/to/working_dir"
RESULT_DIR="path/to/results_dir"
```

The script also takes variables for replication rates *p_repl*, mutation rates *p_mut*, length of sequences *L*, number of sequences *N* and final time of simulation steps *t_final*. To test for different replication and mutation rates, write them in a list:

```
p_repl=(0.5 0.53)
p_mut=(0.0001)
L=100
N=20
t_final=10
```
The script currently runs 20 simulations per parameter set. To change that, replace the value in line:

```
N_SIM=20
```

The script will first run multiple simulations per parameter set and create count files, then create the header files. Afterwards it collects the additional count tables for number of new mutations and number of sequences and puts the average table on top of each simulations collection directory.

The tables for plots are created by **computeNucleotideDiversity.r** in each single simulation directory.

The averaged plots are then created by **averageOverSimulationTables.r**. Averaged plots can be found on top of each simulations collection directory. 

After execution the folder will look for example like this:
```
   .
    ├── ...
    ├── simulations                                                 # Main simulations folder
    │   ├── python_out_p_repl_0.5_p_mut_0.0001                      # Python output tables
    │   ├── python_out_p_repl_0.53_p_mut_0.0001
    │   ├── simulation_p_repl_0.5_p_mut_0.0001                      # Simulation trajectories folder
    │   └── simulation_p_repl_0.53_p_mut_0.0001
    │       ├── python_out_mean_number_of_new_mut_table.txt         # Average of Python tables 
    │       ├── .pdfs                                               # Different averaged plots                       
    │       ├── timesteps_p_repl_0.5_p_mut_0.0001_L_100_N_20_sim_1  # Individual simulation folder
    │       ├── timesteps_p_repl_0.5_p_mut_0.0001_L_100_N_20_sim_2  # Individual simulation folder
    │       ...
    └── ...
```
