# 自动编译和运行
import os, re
import argparse
import json
from pathlib import Path
import os.path as osp
import difflib

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fortran_log", type=str, help="fortran ouput log path")
    parser.add_argument("-c", "--cpp_log", type=str, help="cpp output log path")
    parser.add_argument("-o", "--diff_output", type=str, default="diff.html", help="ouput html file")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
        
    
    with open(args.fortran_log) as f, open(args.cpp_log) as c:
        fout = f.readlines()
        cout = c.readlines()
        diff = difflib.HtmlDiff().make_file(fout, cout, fromdesc="fortran", todesc="c++")
        with open(args.diff_output, 'w') as f:
            f.write(diff)

    
    
    