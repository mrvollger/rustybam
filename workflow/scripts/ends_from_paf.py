#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import os
import sys
import argparse
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("infile", help="positional input")
    parser.add_argument("-d", help="store args.d as true if -d",
        action="store_true", default=False)
    parser.add_argument("-w","--width",help="end size", type=int, default=5e4)
    args = parser.parse_args()
    df = pd.read_csv(args.infile, sep="\t", header=0, usecols=range(12), names=list(range(12)))
    
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


