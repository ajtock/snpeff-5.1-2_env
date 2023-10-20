#!/bin/bash

# Created: 09/05/2023

# Create a self-contained conda environment for
# installation of snpeff

mamba env create --name snpeff \
                 --file snpeff.yaml
