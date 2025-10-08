<div align="center">
    <h1>Distributed Stencil Performance Analysis</h1>
    <h3>Author: Christian Faccio</h3>
    <h5>Email: christianfaccio@outlook.it</h4>
    <h5>Github: <a href="https://github.com/christianfaccio" target="_blank">christianfaccio</a></h5>
    <h6>This project aims at analyzing the performance of distributed stencil computations in a high-performance computing environment. Starting from a serial code implementation of the problem, a mixed approach using MPI and OpenMP is employed to evaluate the benefits of parallelism.</h6>
</div>

---

<div align="center">
   <img src="slides/assets/energy_evolution_periodic.gif" alt="Energy Evolution">
</div>

## Description

The objective of this project is to evaluate the performance of a parallel implementation of a 2D stencil computation using a hybrid MPI and OpenMP approach. Usually used in numerical simulations, stencil computations involve updating each point in a grid based on the values of its neighboring points. Here we consider the **heat equation** as a case study:

```math
\frac{\partial u(t,\vec{x})}{\partial t} = \alpha \nabla^2 u(t,\vec{x})
```

where $\alpha$ is the heat diffusivity and $u$ the energy at position $\vec{x}$ and time $t$. Discretizing the equation in two dimensions  on an $m\times l$ grid using finite differences leads to the following update rule for each grid point:

```math
U_{m,l}^{n+1} = U_{m,l}^n + \frac{\alpha \Delta t}{\Delta x^2} \left( U_{m-1,l}^n + U_{m+1,l}^n - 2U_{m,l}^n \right) + \frac{\alpha \Delta t}{\Delta y^2} \left( U_{m,l-1}^n + U_{m, l+1}^n - 2U_{m,l}^n \right) 
```

This approximation is known as the **five-points stencil**. The project involves implementing both a serial and a parallel version of this computation, analyzing their performance in terms of execution time and scalability.

## Project Structure

```
.
├── assignment.pdf          # Project assignment
├── include                 # Header files
│   ├── parallel.h
│   └── serial.h
├── Makefile                # Makefile to build the project
├── output                  # Output files (.out from Leonardo)
│   ├── strong_scaling
│   ├── threads_scaling
│   └── weak_scaling
├── plots                   # Plots generated from the output data  
├── python                  # Python scripts for data analysis and visualization  
│   ├── parse_results.py
│   ├── plot_results.py
│   └── visualize_grid.py
├── README.md               # Project documentation
├── results                 # Results of the experiments (.csv files)
├── requirements.txt        # Project dependencies
├── scripts                 # Scripts for running experiments on Leaonardo
└── src                     # Source code files
    ├── parallel.c
    └── serial.c
```

## How to build and run

First of all, clone the repository:

```bash
git clone git@github.com:christianfaccio/stencil_method.git
```

(eventually) create and activate a virtual environment:

```bash
uv venv stencil
source stencil/bin/activate
uv pip install -r requirements.txt
```

Then, build and run the project using the provided Makefile:
```
make serial <args>   # For the serial version
make parallel <args> # For the parallel version
```

To clean the build files, use:
```
make clean_all
```

## Leonardo Cluster

This project is designed to be used on a high-performance computing cluster, such that an analysis of the performance of the parallel implementation can be conducted. The cluster used for this project is Leonardo, specifically the DCGP partition.

Each DCGP node features:

- 2x Intel(R) Xeon(R) Platinum 8480+ CPUs
- 112 physical cores in total (56 per socket)
- Hyper-Threading disabled (1 thread per core)
- A complex NUMA architecture with 8 nodes.

For everyday login, use the following command:

```
step ssh login '<email>' --provisioner cineca-hpc # login

ssh <username>@login.leonardo.cineca.it
```

if a warning (`WARNING: REMOTE HOST IDENTIFICATION HAS CHANGED!`) appears, just run the following command and login:
```
ssh-keygen -f ~/.ssh/known_hosts -R login.leonardo.cineca.it; for keyal in ssh-rsa ecdsa-sha2-nistp256; do for address in login01-ext.leonardo.cineca.it login02-ext.leonardo.cineca.it login05-ext.leonardo.cineca.it login07-ext.leonardo.cineca.it; do ssh-keyscan -t  ${keyal} ${address} | sed "s/\b${address}/login*.leonardo.cineca.it/g" >> ~/.ssh/known_hosts; done; done
ssh <username>@login.leonardo.cineca.it
```

`saldo -b --dcgp` tells how many core-hours (budget) is left for the **group**.

To share file TO/FROM local machine and Leonardo, I used the `scp` command via datamovers:

```
scp path/to/local/file <username>@data.leonardo.cineca.it:/leonardo/home/<username>
scp <username>@data.leonardo.cineca.it:/leonardo/home/<username>/path/to/remote/file path/to/local/
```

Finally, to submit jobs to the DCGP partition, use the `sbatch` command:

```
# Assure no modules
module purge
# Load openmpi module (leonardo)
module load openmpi/4.1.6--gcc--12.2.0
# Compile
scripts/compile.sh

# Submit job
sbatch scripts/threads_scaling.sh
sbatch scripts/strong_scaling.sh
sbatch scripts/weak_scaling.sh
```

> [!WARNING]
> The `go_dcgp.sbatch` script is used internally by the other .sh scripts so don't use it directly (at least for this project).

## Choices made

This section is made to explain why certain code chunks were implemented in a specific way or the idea behind some configurations.

- **Separate context**: `MPI_Comm_dup(MPI_COMM_WORLD, &STENCIL_WORLD)` is best practice to isolate the parallel code from other tasks in the computer;
- `MPI_THREAD_FUNNELED` is used because only the main thread will make MPI calls;
- **Domain Decomposition**: Total domain is divided among processes, using an heuristics to avoid creating thin slices (with form factor) and using a simple factorization algorithm.
- **Sources init**: Task 0 assigns randomly the sources to the tasks, shares with them using `MPI_Bcast` and then each task randomly places the sources in its subdomain.
- **MPI Communication**: first (SEND) part of the buffer is filled with local data, then it is shared and saved in the second part of the buffer (RECV) such that all the tasks know the correct values to use in the stencil computation.
- **Asynchronous communication**: `MPI_Isend` and `MPI_Irecv` are used to overlap communication and computation, improving performance.
- **OpenMP Parallel Regions**: `#pragma GCC unroll 4` is used to suggest the compiler to unroll the loop 4 times, potentially improving performance by reducing loop overhead and increasing instruction-level parallelism. `schedule(static)` is used to divide the loop iterations into contiguous chunks, which can improve cache performance and reduce overhead. `export OMP_PLACES=cores` and `export OMP_PROC_BIND=close` are set to bind threads to specific cores, reducing context switching and improving cache locality. Finally, `MPI_Reduce` is used to compute the total energy. First touch policy is used to allocate memory in a NUMA-aware manner, improving memory access times.
- **Optimizations**: `register` keyword is used to suggest the compiler to store variables in CPU registers for faster access. `-O3` optimization flag is used to enable high-level optimizations during compilation, potentially improving performance. `-march=native` enables the use of all instruction sets available on the host machine, optimizing the code for the specific architecture of the CPU.
- **Why 14 threads per task?**: Chosing 4 would be better in terms of cache hits, but on a HPC cluster where bigger problems are run, this would lead to poorer results. 8 threads per task would have lead to a strange division of tasks per node and so 14 has been chosen as a compromise.