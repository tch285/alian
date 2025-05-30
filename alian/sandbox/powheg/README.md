# how to use powheg + pythia

- some instructions at https://github.com/matplo/heppyy/blob/main/heppyy/example/powheg_pythia/README.md 
- here we added Zjet config (note: for muonic decays) additional Zjet that requires 
- this relies on heppyy installed (https://github.com/matplo/heppyy) as the rest of the alian

# two-steps:

1. generate .lhe file with powheg
2. shower+hadronize (analyze) with pythia8
```
if [ -e "pwgevents.lhe" ]; then
	$HEPPYY_DEV/heppyy/example/powheg_pythia/pythia8_powheg.py pwgevents.lhe --nev 1000 --verbose
fi
```

## one liner

- to run Zjet look into `run_Zj.sh` or just run it... similar example for b's `run_b.sh`

# useful code example

- load lhe .lhe file to pythia and run for some events (see [pythia8_powheg.py](https://github.com/matplo/heppyy/blob/main/heppyy/example/powheg_pythia/pythia8_powheg.py))
- we read .lhe file + run pythia and fastjet for each event... 

```
./pythia8_powheg.py events.lhe --nev 1000 --verbose
```

