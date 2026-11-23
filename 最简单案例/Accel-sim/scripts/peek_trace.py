#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Show Accel-Sim kernelslist.g and the first kernel header + a few insts."""
from __future__ import print_function
import os
import sys


def main():
    if len(sys.argv) < 2:
        print("usage: python peek_trace.py <kernelslist.g>")
        sys.exit(1)
    lst = sys.argv[1]
    folder = os.path.dirname(os.path.abspath(lst))
    print("kernelslist:", lst)
    print("-" * 60)
    for line in open(lst, "r"):
        line = line.strip()
        if not line:
            continue
        print("  ", line)
        if line.startswith("kernel-") or line.endswith(".traceg") or line.endswith(".tracez"):
            kpath = os.path.join(folder, line)
            if os.path.isfile(kpath):
                print("      ->", kpath)
                n = 0
                for kl in open(kpath, "r"):
                    print("      ", kl.rstrip())
                    n += 1
                    if n >= 24:
                        print("      ...")
                        break
    print()
    print("This is a SASS instruction trace (what the GPU executed).")
    print("It is NOT gem5 stats. Accel-Sim replays it, AccelWattch then emits Watts.")


if __name__ == "__main__":
    main()
