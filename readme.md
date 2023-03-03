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
- ModelingToolkit
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

