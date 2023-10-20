#!/bin/bash

# Created: 09/05/2023
# Updated: 20/09/2023 (bcftools=1.17)

# Update self-contained conda environment for
# installation of snpeff

mamba env update --name snpeff \
                 --file snpeff.yaml
