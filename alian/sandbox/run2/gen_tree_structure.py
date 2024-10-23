#!/usr/bin/env python

import ROOT
import yaml
import argparse

def main():
	parser = argparse.ArgumentParser(description='read a root file and prepare a yaml file with tree structure')
	parser.add_argument('input', help='Input file')
	parser.add_argument('-o', '--output', help='Output file', default='tree_structure.yaml')
 
	args = parser.parse_args()

	input_file = ROOT.TFile(args.input)
	
	# find all the keys in the file - recurse into TDrirectoryFile if necessary
	# find all the trees
	# build a dictionary with the tree structure
	# write the dictionary to a yaml file
 
	keys = input_file.GetListOfKeys()
	tree_dict = {}
	for key in keys:
		obj = key.ReadObj()
		if isinstance(obj, ROOT.TTree):
			tree_dict[obj.GetName()] = {}
			tree = obj
			branches = tree.GetListOfBranches()
			for branch in branches:
				branch_name = branch.GetName()
				leaf = branch.GetLeaf(branch_name)
				leaf_type = leaf.GetTypeName()
				tree_dict[obj.GetName()][branch_name] = leaf_type
		elif isinstance(obj, ROOT.TDirectoryFile):
			dir_dict = {}
			for key in obj.GetListOfKeys():
				dir_obj = key.ReadObj()
				if isinstance(dir_obj, ROOT.TTree):
					dir_dict[dir_obj.GetName()] = {}
					tree = dir_obj
					branches = tree.GetListOfBranches()
					for branch in branches:
						branch_name = branch.GetName()
						leaf = branch.GetLeaf(branch_name)
						leaf_type = leaf.GetTypeName()
						dir_dict[dir_obj.GetName()][branch_name] = leaf_type
			if len(dir_dict) > 0:
				tree_dict[obj.GetName()] = dir_dict
		else:
			print('Unknown object type: ', obj)
   
   # write the dictionary to a yaml file
	with open(args.output, 'w') as yaml_file:
		yaml.dump(tree_dict, yaml_file, default_flow_style=False)
 
if __name__ == '__main__':
	main()