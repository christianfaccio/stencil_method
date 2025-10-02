<div align="center">
    <h1>Distributed Stencil Performance Analysis</h1>
    <h3>Author: Christian Faccio</h3>
    <h5>Email: christianfaccio@outlook.it</h4>
    <h5>Github: <a href="https://github.com/christianfaccio" target="_blank">christianfaccio</a></h5>
    <h6>This project aims at analyzing the performance of distributed stencil computations in a high-performance computing environment. Starting from a serial code implementation of the problem, a mixed approach using MPI and OpenMP is employed to evaluate the benefits of parallelism.</h6>
</div>

---

## Decisions for parallel code

- #### Use simple for loop with sendrecv instead of multiple groups/communicators for the **gather** operation between processes.
    1. Overhead: Creating groups and communicators has significant overhead, and this would be repeated every iteration
    2. Unnecessary Complexity: The document emphasizes that you should "create a separated context for the application" but only for logical isolation, not for simple point-to-point patterns
    3. Performance: Simple point-to-point communication (especially with MPI_Sendrecv or non-blocking calls) is highly optimized and exactly designed for this halo exchange pattern
    4. Standard Pattern: Your use case is a classic halo exchange or ghost cell exchange pattern, which is one of the most common communication patterns in HPC and is best handled with point-to-point communication

- #### Better to use a **mixed approach** (MPI + OpenMP) instead of having many virtual processes with MPI only. 
    1. Memory Efficiency: Each MPI process has its own memory space, which can lead to high memory overhead when you have many processes. OpenMP threads share the same memory space, which is more memory efficient.
    2. Reduced Communication Overhead: MPI processes communicate via message passing, which can introduce significant overhead, especially with many processes. OpenMP threads can communicate via shared memory, which is generally faster.
    3. Scalability: A mixed approach allows you to scale more effectively on modern multi-core and many-core architectures. You can use MPI to distribute work across nodes and OpenMP to utilize multiple cores within each node.
    4. Load Balancing: Using OpenMP within each MPI process allows for dynamic load balancing among threads, which can be beneficial for irregular workloads.
    5. Simplicity in Code: A mixed approach can simplify the code structure by reducing the number of MPI processes and leveraging OpenMP for parallelism within each process.

- #### For scaling on Leonardo:
    1. **OpenMP** scaling -> 1 MPI process and scale with threads. Find the best thread/task ratio
    2. **MPI** scaling:
        
        1. Strong scaling -> Keep the problem size fixed and increase the number of MPI processes
        2. Weak scaling -> Increase the problem size proportionally to the number of MPI processes






## Explanations

- #### Domain Decomposition (lines 271-297 of src/parallel.c):
  - Calculates a form factor (aspect ratio) of the domain S
  - Decides whether to use 1D or 2D decomposition based on the number of processes
  - For 1D: all processes along the longer dimension
  - For 2D: factorizes Ntasks and distributes processes to match domain aspect ratio
  - Maps rank to 2D grid coordinates using row-major ordering
  - Determines NORTH, SOUTH, EAST, WEST neighbors
  - Handles both periodic and non-periodic boundaries
  - Each process gets base size: s = S[dim] / Grid[dim]
  - Remainder r = S[dim] % Grid[dim]
  - First r processes in each dimension get +1 extra cell: mysize[dim] = s + (coord < r)