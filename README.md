#A general program to perform simultaneous maximum likelihood fittings and simulation of the decay curves with gates on neutron. The decay curves can be given as a form of unbin or bin data.
### To run benchmark:

```bash
cd build-pgenfit
./runbenchmark.sh
```

ignore errors

comparing with mlhfitv1

```bash
cd pgenfit/dev/mlhfitv1_forbenchmark
./runbenchmark.sh
```

- To run simulation

```bash
./simulationC parmsex.txt simparmsex.txt output.root
```

- See [https://nbviewer.org/github/vihophong/pgenfit/blob/pythonize/PSDetSim_230726.ipynb](https://nbviewer.org/github/vihophong/pgenfit/blob/pythonize/PSDetSim_230726.ipynb) some examples of running simulations.
