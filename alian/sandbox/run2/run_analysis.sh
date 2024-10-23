./analysis.py /Users/ploskon/data/alice/run2/file1.root jet_tree_min.yaml -o analysis_results_1.root 
./analysis.py /Users/ploskon/data/alice/run2/file2.root jet_tree_min.yaml -o analysis_results_2.root 
./analysis.py /Users/ploskon/data/alice/run2/file3.root jet_tree_min.yaml -o analysis_results_3.root 
./analysis.py /Users/ploskon/data/alice/run2/file4.root jet_tree_min.yaml -o analysis_results_4.root 
./analysis.py /Users/ploskon/data/alice/run2/file5.root jet_tree_min.yaml -o analysis_results_5.root 

hadd -f analysis_results.root analysis_results_1.root analysis_results_2.root analysis_results_3.root analysis_results_4.root analysis_results_5.root
