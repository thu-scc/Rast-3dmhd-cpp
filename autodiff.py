# 自动编译和运行
import os, re
import argparse
import json
from pathlib import Path
import os.path as osp
import subprocess
from subprocess import Popen, PIPE
from threading import Timer
import difflib
import sys

MPI = "openmpi/dykxdir"
TIME_OUT = 5
F_DIR = None
C_DIR = None
compile_fortran_cmd = [
    f". /opt/spack_new/share/spack/setup-env.sh",
    "spack load intel-oneapi-compilers",
    f"spack load intel-oneapi-mpi/z3z7huc",
    "mpiifort -O0 -g -std=legacy -ffpe-summary='none' -c  3dmhd.f -o 3dmhd.o",
    "mpiifort -O0 -g -std=legacy -ffpe-summary='none' -c  3dmhdset.f -o 3dmhdset.o",
    "mpiifort -O0 -g -std=legacy -ffpe-summary='none' -c  3dmhdsub.f -o 3dmhdsub.o",
    "mpiifort -O0 -g -std=legacy -ffpe-summary='none' 3dmhd.o 3dmhdset.o 3dmhdsub.o   -o 3dmhd"
]
compile_C_cmd = [
    f". /opt/spack_new/share/spack/setup-env.sh",
    f"spack load {MPI}",
    f"mpic++ -O0 -std=c++20 -c 3dmhd.cpp -o 3dmhd.o",
    f"mpic++ -O0 -std=c++20 -o 3dmhd 3dmhd.o"
]



def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fortran_dir", type=str, default="/home/Siunaus/3dmhd/3dmhd", help="image dir to extract images")
    parser.add_argument("-c", "--cpp_dir", type=str, default="/home/Siunaus/3dmhd/3dmhd-cpp", help="prompt dir")
    parser.add_argument("-o", "--diff_output", type=str, default="test.xls", help="输出的xls文件")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    F_DIR = args.fortran_dir
    C_DIR = args.cpp_dir
    # 编译
    if osp.exists(F_DIR) == False or osp.exists(C_DIR) == False:
        print("文件夹不存在")
        exit(0)
    

    subprocess.run("\n".join(compile_fortran_cmd), shell=True, 
                   stdout=subprocess.PIPE,
                   cwd=F_DIR).check_returncode()
    subprocess.run("\n".join(compile_C_cmd), shell=True, 
                   stdout=subprocess.PIPE,
                   cwd=C_DIR).check_returncode()
    print("compile done")
    
    run_fortran_cmd = [
        f". /opt/spack_new/share/spack/setup-env.sh",
        f"spack load intel-oneapi-mpi/z3z7huc",
        "mpirun -np 24 ./3dmhd"
    ]
    run_c_cmd = [
        f". /opt/spack_new/share/spack/setup-env.sh",
        f"spack load {MPI}",
        "mpirun -np 24 ./3dmhd"
    ]
    pf = Popen("\n".join(run_fortran_cmd), stdout=PIPE, stderr=PIPE, shell=True, cwd=F_DIR)
    timerf = Timer(TIME_OUT, pf.kill)
    pc = Popen("\n".join(run_c_cmd), stdout=PIPE, stderr=PIPE, shell=True, cwd=C_DIR)
    timerc = Timer(TIME_OUT, pc.kill)

    try:
        timerf.start()
        timerc.start()
        pf.wait()
        pc.wait()
    except Exception as e:
        print(e)
        timerf.cancel()
        timerc.cancel()
    
    print("run done")
        
    with open("fout.log", "w") as f, open("cout.log", "w") as c:
        f.write(pf.stdout.read().decode())
        c.write(pc.stdout.read().decode())
    
    with open("fout.log") as f, open('cout.log') as c:
        fout = f.readlines()
        cout = c.readlines()
        diff = difflib.HtmlDiff().make_file(fout, cout, fromdesc="fortran", todesc="c++")
        with open('diff.html', 'w') as f:
            f.write(diff)

    
    
    