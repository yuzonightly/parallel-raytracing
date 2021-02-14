# raytracing-discovery

Goal: use PThreads and MPI to parallelize a simple raytracing algorithm.

## Execution

### hybrid_version

```console
me@mars:~/hybrid_version$ mkdir build && cd build
me@mars:~/hybrid_version/build$ cmake ..
me@mars:~/hybrid_version/build$ make
me@mars:~/hybrid_version/build$ mpirun -np [number of processes] ./ray_tracing_one_week --threads [number of threads] --scene earth earth.ppm
```

### sequential_version

```console
me@mars:~/hybrid_version$ mkdir build && cd build
me@mars:~/hybrid_version/build$ cmake ..
me@mars:~/hybrid_version/build$ make
me@mars:~/hybrid_version/build$ ./ray_tracing_one_week --scene earth earth.ppm
```

### mpi_version

```console
me@mars:~/hybrid_version$ mkdir build && cd build
me@mars:~/hybrid_version/build$ cmake ..
me@mars:~/hybrid_version/build$ make
me@mars:~/hybrid_version/build$ mpirun -np [number of processes] ./ray_tracing_one_week --scene earth earth.ppm
```

### thread_version

```console
me@mars:~/hybrid_version$ mkdir build && cd build
me@mars:~/hybrid_version/build$ cmake ..
me@mars:~/hybrid_version/build$ make
me@mars:~/hybrid_version/build$ ./ray_tracing_one_week --threads [number of threads] --scene earth earth.ppm
```
