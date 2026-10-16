#!/usr/bin/env python3
"""Run a CHARMM-GUI/OpenMM system: equilibration, then production in restartable segments.

The system directory is the `openmm` folder produced by CHARMM-GUI (it holds openmm_run.py,
toppar.str, step3_input.psf/.crd, sysinfo.dat, step4_equilibration.inp and step5_production.inp).
It is copied to the replicate directory; the folders that toppar.str refers to (../gdp, ../gtp,
../lig, ../toppar, ...) are copied next to the replicate directory.

Each production segment writes step5_<i>.rst/.dcd and, when it finishes without error,
step5_<i>.done. Segments that already have a marker (and their files) are skipped, and the
partial output of an interrupted segment is deleted before it is rerun, so the script can be
started again after a Colab session ends. The length of a segment is nstep * dt from
step5_production.inp (1 ns for 250,000 steps of 4 fs) and is printed with the progress.

Usage:
    python run_md.py --system-dir <.../charmm-gui-XXXX/openmm> --rep-dir <.../reps_s12/rep1> --n-cycles 50
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import time

PAR = ["-t", "toppar.str", "-p", "step3_input.psf", "-c", "step3_input.crd"]


def read_inp(path):
    """Parse 'key = value  # comment' lines of a CHARMM-GUI/OpenMM .inp file."""
    kv = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("#")[0].strip()
            if "=" in line:
                key, val = line.split("=", 1)
                kv[key.strip()] = val.strip()
    return kv


def segment_ns(inp_path):
    """Length of one production segment in ns (nstep * dt, dt in ps)."""
    kv = read_inp(inp_path)
    return int(kv["nstep"]) * float(kv["dt"]) / 1000.0


def support_dirs(toppar_text):
    """Folders referenced as ../<name>/ in toppar.str."""
    return sorted(set(re.findall(r"\.\./([\w\-\.]+)/", toppar_text)))


def is_done(rep_dir, tag):
    """True if the segment has its marker and its restart and trajectory files."""
    marker = os.path.join(rep_dir, tag + ".done")
    needed = [os.path.join(rep_dir, tag + ext) for ext in (".rst", ".dcd")]
    if not os.path.exists(marker):
        return False
    if not all(os.path.exists(f) for f in needed):
        print(f"[{tag}] marker present but output files missing - rerunning", flush=True)
        return False
    return True


def purge(rep_dir, tag):
    """Delete the output of a segment that did not finish."""
    for ext in (".rst", ".dcd", ".out", ".done"):
        path = os.path.join(rep_dir, tag + ext)
        if os.path.exists(path):
            os.remove(path)
            print(f"[removed incomplete {tag}{ext}]", end=" ", flush=True)


def prepare_rep_dir(system_dir, rep_dir):
    system_dir, rep_dir = os.path.abspath(system_dir), os.path.abspath(rep_dir)
    root, parent = os.path.dirname(system_dir), os.path.dirname(rep_dir)
    os.makedirs(parent, exist_ok=True)
    with open(os.path.join(system_dir, "toppar.str")) as fh:
        wanted = support_dirs(fh.read())
    for name in wanted:
        src, dst = os.path.join(root, name), os.path.join(parent, name)
        if not os.path.isdir(src):
            raise FileNotFoundError(f"toppar.str refers to ../{name}/ but {src} does not exist")
        if not os.path.isdir(dst):
            shutil.copytree(src, dst)
            print(f"copied {name}/ next to the replicate directory", flush=True)
    if not os.path.exists(os.path.join(rep_dir, "openmm_run.py")):
        if os.path.isdir(rep_dir) and os.listdir(rep_dir):
            # never delete a directory that may hold trajectories
            raise RuntimeError(f"{rep_dir} exists but contains no openmm_run.py (incomplete copy?); "
                               "move or remove it by hand and run again")
        if os.path.isdir(rep_dir):
            os.rmdir(rep_dir)
        shutil.copytree(system_dir, rep_dir)
        print(f"created {rep_dir}", flush=True)
    return rep_dir


def run(cmd, cwd, log):
    with open(os.path.join(cwd, log), "w") as fh:
        return subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT).returncode


def environment_report():
    """Print the OpenMM version and the platforms it offers (recorded in the run log)."""
    try:
        import openmm
    except ImportError:
        print("OpenMM is not importable from this interpreter", flush=True)
        return
    platforms = [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]
    print(f"OpenMM {openmm.__version__}; platforms: {', '.join(platforms)}", flush=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--system-dir", required=True, help="CHARMM-GUI openmm folder")
    p.add_argument("--rep-dir", required=True, help="replicate directory (created if missing)")
    p.add_argument("--n-cycles", type=int, default=50, help="number of production segments")
    p.add_argument("--python", default=sys.executable, help="interpreter used to run openmm_run.py")
    a = p.parse_args(argv)

    environment_report()
    rep = prepare_rep_dir(a.system_dir, a.rep_dir)
    seg = segment_ns(os.path.join(rep, "step5_production.inp"))
    base = [a.python, "-u", "openmm_run.py"]

    if not is_done(rep, "step4_equilibration"):
        purge(rep, "step4_equilibration")
        print("equilibration ...", end=" ", flush=True)
        t0 = time.time()
        rc = run(base + ["-i", "step4_equilibration.inp", *PAR, "-b", "sysinfo.dat",
                         "-orst", "step4_equilibration.rst", "-odcd", "step4_equilibration.dcd"],
                 rep, "step4_equilibration.out")
        if rc:
            print("FAILED - see step4_equilibration.out", flush=True)
            return rc
        open(os.path.join(rep, "step4_equilibration.done"), "w").write("ok")
        print(f"done ({(time.time() - t0) / 60:.0f} min)", flush=True)
    else:
        print("equilibration complete - skipped", flush=True)

    for i in range(1, a.n_cycles + 1):
        tag = f"step5_{i}"
        if is_done(rep, tag):
            continue
        purge(rep, tag)
        prev = "step4_equilibration" if i == 1 else f"step5_{i - 1}"
        print(f"production segment {i}/{a.n_cycles} (ends at {i * seg:g} ns) ...", end=" ", flush=True)
        t0 = time.time()
        rc = run(base + ["-i", "step5_production.inp", *PAR, "-irst", f"{prev}.rst", "-orst", f"{tag}.rst",
                         "-odcd", f"{tag}.dcd"], rep, f"{tag}.out")
        if rc:
            print(f"FAILED - see {tag}.out", flush=True)
            return rc
        open(os.path.join(rep, tag + ".done"), "w").write("ok")
        print(f"done ({(time.time() - t0) / 60:.0f} min)", flush=True)

    n = sum(1 for i in range(1, a.n_cycles + 1) if os.path.exists(os.path.join(rep, f"step5_{i}.done")))
    print(f"completed: {n}/{a.n_cycles} segments = {n * seg:g} ns", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
