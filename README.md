SigrefMC experiments
====================
This repository hosts the source code of the experimental section in the paper on [SigrefMC](https://github.com/utwente-fmt/sigrefmc) submitted to TACAS 2016.

You can contact the main author of SigrefMC at <t.vandijk@utwente.nl>.

SigrefMC source code: https://github.com/utwente-fmt/sigrefmc  
Sylvan source code: https://github.com/utwente-fmt/sylvan  

Files
=====
The `sigref` program is a modified version of [sigref 1.5](http://ira.informatik.uni-freiburg.de/software/sigref/).
It is modified to disable quotient extraction and to print time per iteration.

The `sigref_gmp` and `sigref_floating` programs are versions of sigref used in [Wimmer, Becker, Correctness Issues of Symbolic Bisimulation Computation for Markov Chains (2010)](http://link.springer.com/chapter/10.1007%2F978-3-642-12104-3_22).
We use the `sigref_gmp` tool since it uses a modified version of CUDD to implement MTBDDs with rational numbers as leaves.

The `sigrefmc` tool is the version of the SigrefMC tool used to generate the experimental data.

The `sigrefmc_ht` tool is a version of the SigrefMC tool that uses a hash table instead of a skip list in `refine`.

The `models` directory contains the benchmark models included in the distributions of above tools.

The `out` directory contains all log files generated for the experiments.

The `sigrefmc` directory contains the source code of the tool.

The `recompile.py` script recompiles `sigrefmc`.

The `exp.py` file runs the experiments and reports data.

Experiments
===========
To reproduce the results, start with a clean `out` directory and run:

```
exp.py run
```

To view the results of experiments and generate latex tables written to `results_*.tex`, run:

```
exp.py report
```

The script runs the following commands for LTS strong bisimulation:
```
sigref --bisi=2 --verbosity=1 --infile=<model>
sigrefmc -b 2 -w <workers> <model>
```

The script runs the following commands for LTS branching bisimulation:
```
sigref --bisi=1 --verbosity=1 --infile=<model>
sigrefmc -b 1 -w <workers> <model>
```

The script runs the following commands for CTMC strong bisimulation:
```
sigref_gmp <model>
sigrefmc -l fraction -w <workers> <model>

```
The script runs the following commands for IMC strong bisimulation:
```
sigref --bisi=2 --verbosity=1 --infile=<model>
sigrefmc -b 2 -w <workers> <model>
```

The script runs the following commands for IMC branching bisimulation:
```
sigref --bisi=1 --verbosity=1 --infile=<model>
sigrefmc -l floating -b 1 -w <workers> <model>
```
