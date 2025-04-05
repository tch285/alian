#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This script generates a pwgseeds.dat file with unique random seeds.
# It is intended for use with POWHEG

import random

def generate_pwgseeds(filename="pwgseeds.dat", num_seeds=10):
    """
    Generate a pwgseeds.dat file with unique random seeds.

    :param filename: Name of the output file (default: pwgseeds.dat)
    :param num_seeds: Number of unique random seeds to generate (default: 10)
    """
    seeds = set()
    while len(seeds) < num_seeds:
        seed = random.randint(1, 999999)
        seeds.add(seed)

    with open(filename, "w") as f:
        for seed in seeds:
            f.write(f"{seed}\n")
    print(f"{filename} generated with {num_seeds} unique random seeds.")

# Customize the number of seeds if needed
if __name__ == "__main__":
    generate_pwgseeds(num_seeds=100000)