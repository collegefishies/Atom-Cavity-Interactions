# Atom-Cavity-Interactions

This is a library for simulating atom-cavity interactions using Python.

## Dependencies

- qutip
- numpy
- pandas

For the Julia, script
- julia
- pyJulia
- QuantumCumulants
- ModelingToolkit@v8.21.0
- OrdinaryDiffEq
- PyPlot
- Tricks


## Files

1. **interaction.py**: This is the library that contains the simulator code.
2. **interaction_picture_cli.py**: This wraps that library in a CLI interface.
3. **interaction_picture_batch_semaphore.py**: This allows for parallelization of the cli picture for sweeping through parameters.
4. **cavity_squeezing.jl**: This is a Julia Script written by Christoph Hotter. It can run much faster code using the Heisenberg picture and their QuantumCumulants. Julia has a massive speedup, if you keep the interpreter running it open. Running the script fresh takes a lot of time.
## Usage

To run the code, navigate to the project directory and use the following command:

```bash
python interaction_picture_batch_semaphore.py -N 4 -d 'np.linspace(-1,1,10)' -g 33 --kappa 1 --gamma 0.36 --omega 10 --driving_strength 0.1 --scale 1 -T 60 -o batch/data --spin_state_command 'qutip.tensor([u,d,u,d])' -p 1
```

This command will generate a batch of .csv files and compress them altogether at the end. Please ensure there is only one `.csv`file per folder, as it will erase any other `.csv` files in that folder.

## Output

The output of the example command will be a batch of .csv files. The total time taken to complete the process will also be displayed.

## Christoph's Julia Packages

Run this commands in Julia after typing the bracket `]`
```
add CollectiveSpins@0.1.5
add Colors@0.12.10
add Combinatorics@1.0.2
add DiffEqCallbacks@2.25.0
add DifferentialEquations@7.1.0
add FileIO@1.16.0
add FixedPointNumbers@0.8.4
add IJulia@1.24.0
add ImageMagick@1.2.2
add Interpolations@0.14.7
add JLD2@0.4.30
add Juno@0.8.4
add Latexify@0.15.18
add LsqFit@0.13.0
add MacroTools@0.5.10
add ModelingToolkit@8.21.0
add OrdinaryDiffEq@6.11.2
add Plots@1.38.5
add ProgressMeter@1.7.2
add PyCall@1.95.1
add PyPlot@2.11.0
add QuantumCumulants@0.2.13
add QuantumOptics@1.0.8
add QuantumOpticsBase@0.3.8
add SteadyStateDiffEq@1.12.0
add StochasticDiffEq@6.44.0
add SymbolicUtils@0.19.11
add Symbolics@4.14.0
add WignerSymbols@2.0.0```
