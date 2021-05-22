#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import os
import sys
import argparse
import pandas as pd
"""
Col	Type	Description
1	string	Query sequence name
2	int	Query sequence length
3	int	Query start (0-based; BED-like; closed)
4	int	Query end (0-based; BED-like; open)
5	char	Relative strand: "+" or "-"
6	string	Target sequence name
7	int	Target sequence length
8	int	Target start on original strand (0-based)
9	int	Target end on original strand (0-based)
10	int	Number of residue matches
11	int	Alignment block length
12	int	Mapping quality (0-255; 255 for missing)
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="positional input")
    parser.add_argument("-d", help="store args.d as true if -d",
        action="store_true", default=False)
    parser.add_argument("--minwidth",help="minimum alignment size to keep", type=int, default=5e4)
    parser.add_argument("-w","--width",help="end size", type=int, default=1000)
    args = parser.parse_args()
    df = pd.read_csv(args.infile, sep="\t", header=0, usecols=range(12), names=list(range(12)))
    
    # filter alignments that are too small to consider  
    # unless they make up the entire contig more or less
    #remove = (df[10] < args.minwidth) & (df[1] > 2*args.minwidth)
    #df = df[~remove]
    
    outfmt="{}\t{}\t{}\t{}\t{}"
    for name, g in df.groupby(0):
        length = g[1].min()
        mmin = g[2].min()
        mmax = g[3].max()
        if (mmax-mmin) < 2*args.width:
            print(outfmt.format(name, mmin, mmax, length,"full"))
        else:
            print(outfmt.format(name, mmin, args.width,length, "start"))
            print(outfmt.format(name, mmax-args.width, mmax, length, "end"))


