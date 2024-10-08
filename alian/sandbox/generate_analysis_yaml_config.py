#!/usr/bin/env python

# open ROOT file using uproot
# list the trees
# list branches in the trees
# write a yaml file which would describe the structure of the trees

import uproot
import yaml
import argparse

def print_root_keys(file_path):
    file = uproot.open(file_path)
    print("Keys in the ROOT file:")
    for key in file.keys(recursive=True):  # Use recursive=True to see keys in subdirectories
        print('- key:', key, type(key))

# Function to recursively list trees in the ROOT file
def find_trees(root_directory):
    trees = {}
    # Loop over keys in the directory
    for key, item in root_directory.items():
        # print(' - key:', key, type(key), type(item))
        # If the item is a tree, add it to the dictionary
        if isinstance(item, uproot.models.TTree.Model_TTree_v20):
        # if 'uproot.models.TTree.Model_TTree' in str(type(item)):
            print(' - - found tree:', key)
        # Use the key object to get the cycle
            cycle = root_directory.key(key).fCycle
            print(' - - - current cycle:', cycle)
            key_no_cycle = key.split(';')[0]
            if trees.get(key_no_cycle):
                if trees[key_no_cycle][1] < cycle:
                    print(' - - - updating cycle:', trees[key_no_cycle][1], cycle)
                    trees[key_no_cycle] = (item, cycle)
                else:
                    print(' - - - keeping cycle:', trees[key_no_cycle][1], 'instead of', cycle)
            else:
                trees[key_no_cycle] = (item, cycle)
        # If the item is a directory, recursively find trees within it
        elif isinstance(item, uproot.reading.ReadOnlyDirectory):
            print(' - updating trees:', item)
            trees.update(find_trees(item))
    trees = {key: value[0] for key, value in trees.items()}
    return trees

# Function to describe the structure of the ROOT file
def describe_root_file(file_path, tree_name=None):
    # Open the ROOT file
    file = uproot.open(file_path)

    # Get the list of trees in the file
    trees = find_trees(file)

    if not trees:
        raise ValueError("No trees found in the file.")

    # If a specific tree is requested, check if it exists
    if tree_name:
        if tree_name not in trees:
            raise ValueError(f"No tree found with name: {tree_name}")
        trees = {tree_name: trees[tree_name]}  # Keep only the requested tree

    # Structure dictionary to hold the YAML structure
    structure = {}

    # Iterate over trees and list branches
    for tree_name, tree in trees.items():        
        structure[tree_name] = {
            'branches': [branch.name for branch in tree.branches]
        }

    return structure

# Function to write structure to a YAML file
def write_yaml(structure, yaml_file_path):
    with open(yaml_file_path, 'w') as yaml_file:
        yaml.dump(structure, yaml_file, default_flow_style=False)

# Main function using argument parser
def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Describe ROOT file structure and export it to YAML.")

    parser.add_argument(
        "root_file",
        type=str,
        help="Path to the ROOT file."
    )

    parser.add_argument(
        "-t", "--tree",
        type=str,
        help="Optional specific tree name to describe.",
        default=None
    )

    parser.add_argument(
        "-o", "--output",
        type=str,
        help="Path to the output YAML file (default: 'tree_structure.yaml').",
        default="tree_structure.yaml"
    )

    # Parse the arguments
    args = parser.parse_args()

    try:
        print_root_keys(args.root_file)

        # Get the ROOT file structure
        structure = describe_root_file(args.root_file, args.tree)

        # Write the structure to a YAML file
        write_yaml(structure, args.output)

        print(f"Structure saved to {args.output}")
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
