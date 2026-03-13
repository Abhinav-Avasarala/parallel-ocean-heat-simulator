# Parallel Ocean Simulation

A C/OpenMP project that simulates ocean surface temperature updates on a 2D grid using iterative stencil computation.

## Overview
This project explores how different work-partitioning strategies affect performance in a shared-memory parallel program. The simulation updates each interior grid cell using the average of the cell and its north, south, east, and west neighbors from the previous time step.

## Implementations
- `ocean_row.c`: row-based static partitioning
- `ocean_col.c`: column-based static partitioning
- `ocean_sq.c`: square/rectangle-based partitioning
- `ocean_row_bsub`: batch runner for row-based benchmarks
- `ocean_col_bsub`: batch runner for column-based benchmarks
- `ocean_sq_bsub`: batch runner for square/rectangle benchmarks

## Benchmarking
The project measures execution time across multiple thread counts:
- 1 thread
- 2 threads
- 4 threads
- 8 threads
- 16 threads

For each configuration, the benchmark runs 10 times and records:
- raw runtime data
- mean runtime
- speedup relative to the 1-thread baseline

## Build
```bash
gcc -O2 -fopenmp -o ocean_row ocean_row.c
gcc -O2 -fopenmp -o ocean_col ocean_col.c
gcc -O2 -fopenmp -o ocean_sq ocean_sq.c

## Example local run
./ocean_row 1026 5000 grid1026_hotTop.txt output.txt

## Example Cluster run
bsub < ocean_row_bsub
bsub < ocean_col_bsub
bsub < ocean_sq_bsub