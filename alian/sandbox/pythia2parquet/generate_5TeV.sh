# ./pythia2parquet.py --py-hardQCD --py-bias --py-ecm 5020 --nev 200000

# ./pythia2parquet.py --py-hardQCD --py-bias --py-ecm 5020 --nev 200000 --py-biasref 200.

nev=10000
./pythia2parquet.py --py-hardQCD --py-bias --py-ecm 5020 --nev $nev --py-pthatmin 20. -o "any_5TeV_pthat20.parquet"
