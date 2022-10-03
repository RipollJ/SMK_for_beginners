#!/bin/python

from contextlib import redirect_stdout

with open(snakemake.output[0], 'w') as f:
    with redirect_stdout(f):
        print(snakemake.input)
        print(type(snakemake.input))
