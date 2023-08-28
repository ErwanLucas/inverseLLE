# inverseLLE
Genetic-algorithm (GA)-based optimisation of resonator dispersion for Kerr comb state tailoring via the Lugiato-Lefever Equation (LLE). Check the Nature Photonics paper (DOI [10.1038/s41566-023-01252-7](https://doi.org/10.1038/s41566-023-01252-7) or the [arxiv version](https://arxiv.org/abs/2209.10294) for more information.

## Running the optimization
This is done by launching the script `lle_dispersion_genetic_optimize.m`.
Two optimization scenarii are demonstrated:

* Flat comb with a tageted power per line over a given bandwidth. This is currently the default scenario (see `run init_polyCombTarget.m` in the file). The initialization `init_polyCombTarget.m` can be modified to set a precise comb shape target and the fitness can be modified that the GA fits that target.
* Octave-spaced dispersive wave uncomment the line `% run init_octaveDW.m` in `lle_dispersion_genetic_optimize.m` and the opmitization will seek to enhance the comb power in a specified comb line.

## Analysing the optimization results
By default the genetic evolution results are stored in a `.mat` file in the `Data` folder. You can retrieve the optimum spectrum and dispersion profile by running `analyze_genetic_dispersion_results.m`. Additional options can be performed such as simulating a laser scan to check if the state can be reached by scanning the detuning across resonance.
