# Hashira
A GPU accelerated PPP solver written in Julia. <br>
Prepare input and position files (examples provided) and then run
```
julia --threads=(num_threads) (path/to/hashira)/main.jl
```
Requirements: LinearAlgebra, PyPlot (which needs PyCall), PerceptualColourMaps and CUDA (only if GPU acceleration is needed)
