import os
import shutil
import subprocess
import sys
from hatchling.builders.hooks.plugin.interface import BuildHookInterface

COMPILE_SCRIPTS = {"compile.sh", "compileMPI.sh"}

DATA_DIRS = ["docs", "EXMP1", "EXMP2", "EXMP3", "EXMP4",
             "EXMP5", "EXMP6", "EXMP7", "EXMP8", "EXMP9", "EXMPA", "EXMPB"]

class CustomBuildHook(BuildHookInterface):
    def initialize(self, version, build_data):
        root = self.root  # repo root
        src_dir = os.path.join(root, "src")
        bin_dir = os.path.join(root, "bin")
        pkg_dir = os.path.join(root, "pecube")

        # Compile Fortran sources
        os.makedirs(bin_dir, exist_ok=True)
        result = subprocess.run(["make", "all"], cwd=src_dir)
        if result.returncode != 0:
            sys.exit(
                "\n[pecube] Fortran compilation failed.\n"
                "Ensure gfortran and MPI (mpif90, mpif77, mpicc) are on PATH.\n"
                "On HPC: module load gcc openmpi  (or equivalent)\n"
            )


        # Copy binaries into the Python package
        pkg_bin = os.path.join(pkg_dir, "bin")
        os.makedirs(pkg_bin, exist_ok=True)
        for name in os.listdir(bin_dir):
            if name in COMPILE_SCRIPTS:
                continue
            src = os.path.join(bin_dir, name)
            if os.path.isfile(src):
                shutil.copy2(src, pkg_bin)
                os.chmod(os.path.join(pkg_bin, name),
                         os.stat(src).st_mode | 0o111)
        build_data["shared_scripts"] = {
            os.path.join(pkg_bin, name): name
            for name in os.listdir(pkg_bin)
            if not name.endswith('.sh')
        }

        # Copy data directories into the Python package
        for name in DATA_DIRS:
            src = os.path.join(root, name)
            dst = os.path.join(pkg_dir, name)
            if os.path.isdir(src):
                shutil.copytree(src, dst, dirs_exist_ok=True)
