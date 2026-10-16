"""Tests for md/run_md.py (with a mock openmm_run.py) and md/check_system.py (synthetic CHARMM-GUI files).

Run:  python -m pytest tests/   (or python tests/test_md_scripts.py)
"""
import os
import shutil
import sys
import tempfile
import textwrap

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "md"))
import check_system as C  # noqa: E402
import run_md as R  # noqa: E402

FAKE_OPENMM_RUN = textwrap.dedent('''
    import argparse, os, sys
    ap = argparse.ArgumentParser()
    for a in ("-i", "-t", "-p", "-c", "-b", "-irst", "-orst", "-odcd"):
        ap.add_argument(a)
    args = ap.parse_args()
    tag = os.path.basename(args.orst)[:-4]
    with open("calls.log", "a") as fh:
        fh.write(f"{tag} irst={args.irst} inp={args.i}\\n")
    open(args.odcd, "a").write("partial")          # output exists before the step finishes (appended, so leftovers show)
    if os.environ.get("FAIL_AT") == tag:
        sys.exit(3)
    open(args.orst, "w").write("rst")
''')


def make_system(root, nstep=250000, dt=0.004, ligand="LIG", hmr=True, charge_shift=0.0, extra_inp=""):
    omm = os.path.join(root, "openmm")
    os.makedirs(omm)
    os.makedirs(os.path.join(root, "gromacs"))
    os.makedirs(os.path.join(root, "gdp"))
    os.makedirs(os.path.join(root, "toppar"))
    open(os.path.join(root, "gdp", "gdp.str"), "w").write("gdp\n")
    open(os.path.join(root, "toppar", "par.prm"), "w").write("par\n")
    open(os.path.join(omm, "openmm_run.py"), "w").write(FAKE_OPENMM_RUN)
    open(os.path.join(omm, "toppar.str"), "w").write("../toppar/par.prm\n../gdp/gdp.str\n")
    for f in ("step3_input.crd", "sysinfo.dat"):
        open(os.path.join(omm, f), "w").write("x\n")
    open(os.path.join(omm, "step4_equilibration.inp"), "w").write("dt = 0.001  # ps\n")
    open(os.path.join(omm, "step5_production.inp"), "w").write(
        f"nstep = {nstep}  # steps\ndt = {dt}  # ps\ntemp = 310.15\nvdw = Force-switch\ncoulomb = PME\n{extra_inp}")
    hm = 3.024 if hmr else 1.008
    rows = [("PROA", "ALA", "N", -0.30, 14.007), ("PROA", "ALA", "HT1", 0.30, hm), ("PROA", "ALA", "HT2", 0.30, hm),
            ("PROB", "GLY", "HN", 0.31, hm), (ligand, ligand, "C1", -1.0 + charge_shift, 12.011),
            (ligand, ligand, "H1", 0.40, hm), ("SOLV", "TIP3", "OH2", -0.834, 15.9994),
            ("SOLV", "TIP3", "H1", 0.417, 1.008), ("SOLV", "TIP3", "H2", 0.417, 1.008)]
    total = sum(r[3] for r in rows) - charge_shift           # the ion neutralises the unshifted system
    rows.append(("IONS", "SOD", "SOD", -total, 22.99))
    lines = ["PSF EXT CMAP CHEQ XPLOR", "", "         1 !NTITLE", " REMARKS x", "", f"{len(rows):10d} !NATOM"]
    for i, (seg, res, name, q, m) in enumerate(rows, 1):
        lines.append(f"{i:10d} {seg:<8s} {1:<8d} {res:<8s} {name:<8s} X {q:14.6f} {m:13.4f} 0")
    lines += ["", "         1 !NBOND: bonds"]
    open(os.path.join(omm, "step3_input.psf"), "w").write("\n".join(lines) + "\n")
    open(os.path.join(root, "gromacs", "step3_input.gro"), "w").write(
        "t\n%d\n" % len(rows) + " x\n" * len(rows) + "  12.4 8.2 8.2\n")
    return omm, len(rows)


def calls(rep):
    return [l.split() for l in open(os.path.join(rep, "calls.log")).read().splitlines()]


def test_segment_length_and_support_dirs():
    with tempfile.TemporaryDirectory() as d:
        omm, _ = make_system(os.path.join(d, "sys"))
        assert R.segment_ns(os.path.join(omm, "step5_production.inp")) == 1.0
        assert R.support_dirs("../toppar/a.prm\n../gdp/g.str\n../gdp/h.str\n") == ["gdp", "toppar"]


def test_full_run_order_and_restart_chain():
    with tempfile.TemporaryDirectory() as d:
        omm, _ = make_system(os.path.join(d, "sys"))
        rep = os.path.join(d, "reps", "rep1")
        assert R.main(["--system-dir", omm, "--rep-dir", rep, "--n-cycles", "4"]) == 0
        c = calls(rep)
        assert [x[0] for x in c] == ["step4_equilibration", "step5_1", "step5_2", "step5_3", "step5_4"]
        assert c[0][1] == "irst=None"
        assert [x[1] for x in c[1:]] == ["irst=step4_equilibration.rst", "irst=step5_1.rst",
                                        "irst=step5_2.rst", "irst=step5_3.rst"]
        assert c[0][2] == "inp=step4_equilibration.inp" and c[1][2] == "inp=step5_production.inp"
        for tag in ["step4_equilibration"] + [f"step5_{i}" for i in range(1, 5)]:
            assert os.path.exists(os.path.join(rep, tag + ".done")), tag
        assert os.path.exists(os.path.join(d, "reps", "gdp", "gdp.str"))   # support folders copied next to rep
        assert os.path.exists(os.path.join(d, "reps", "toppar", "par.prm"))


def test_interrupted_segment_is_cleaned_and_rerun():
    with tempfile.TemporaryDirectory() as d:
        omm, _ = make_system(os.path.join(d, "sys"))
        rep = os.path.join(d, "reps", "rep1")
        os.environ["FAIL_AT"] = "step5_3"
        try:
            assert R.main(["--system-dir", omm, "--rep-dir", rep, "--n-cycles", "4"]) == 3
        finally:
            del os.environ["FAIL_AT"]
        assert os.path.exists(os.path.join(rep, "step5_2.done"))
        assert not os.path.exists(os.path.join(rep, "step5_3.done"))
        assert os.path.exists(os.path.join(rep, "step5_3.dcd"))              # partial output left behind
        assert R.main(["--system-dir", omm, "--rep-dir", rep, "--n-cycles", "4"]) == 0
        tags = [x[0] for x in calls(rep)]
        assert tags.count("step4_equilibration") == 1                         # not repeated
        assert tags.count("step5_1") == 1 and tags.count("step5_2") == 1      # finished segments not repeated
        assert tags.count("step5_3") == 2 and tags.count("step5_4") == 1      # failed one rerun once
        assert open(os.path.join(rep, "step5_3.rst")).read() == "rst"
        assert open(os.path.join(rep, "step5_3.dcd")).read() == "partial"     # leftover of the failed run was deleted


def test_marker_without_files_is_rerun():
    with tempfile.TemporaryDirectory() as d:
        omm, _ = make_system(os.path.join(d, "sys"))
        rep = os.path.join(d, "reps", "rep1")
        assert R.main(["--system-dir", omm, "--rep-dir", rep, "--n-cycles", "2"]) == 0
        os.remove(os.path.join(rep, "step5_2.dcd"))                           # marker stays
        assert R.main(["--system-dir", omm, "--rep-dir", rep, "--n-cycles", "2"]) == 0
        assert [x[0] for x in calls(rep)].count("step5_2") == 2
        assert os.path.exists(os.path.join(rep, "step5_2.dcd"))


def test_missing_support_folder_is_reported():
    with tempfile.TemporaryDirectory() as d:
        omm, _ = make_system(os.path.join(d, "sys"))
        shutil.rmtree(os.path.join(d, "sys", "gdp"))
        try:
            R.main(["--system-dir", omm, "--rep-dir", os.path.join(d, "reps", "rep1")])
        except FileNotFoundError as e:
            assert "gdp" in str(e)
        else:
            raise AssertionError("expected FileNotFoundError")


def test_incomplete_replicate_directory_is_never_deleted():
    with tempfile.TemporaryDirectory() as d:
        omm, _ = make_system(os.path.join(d, "sys"))
        rep = os.path.join(d, "reps", "rep1")
        os.makedirs(rep)
        keep = os.path.join(rep, "step5_1.dcd")
        open(keep, "w").write("precious trajectory")
        try:
            R.main(["--system-dir", omm, "--rep-dir", rep])
        except RuntimeError as e:
            assert "openmm_run.py" in str(e)
        else:
            raise AssertionError("expected RuntimeError")
        assert open(keep).read() == "precious trajectory"
        os.remove(keep)                                                       # an empty directory is fine
        assert R.main(["--system-dir", omm, "--rep-dir", rep, "--n-cycles", "1"]) == 0


def test_check_system_passes_on_good_system():
    with tempfile.TemporaryDirectory() as d:
        root = os.path.join(d, "sys")
        _, n = make_system(root)
        assert C.check_system(root, expect_atoms=n, expect_box=(124, 82, 82), log=lambda *_: None)


def test_check_system_flags_problems():
    quiet = lambda *_: None
    with tempfile.TemporaryDirectory() as d:
        root = os.path.join(d, "a")
        _, n = make_system(root, hmr=False)
        assert not C.check_system(root, log=quiet)                            # hydrogens not repartitioned
        root = os.path.join(d, "b")
        _, n = make_system(root, charge_shift=0.5)
        assert not C.check_system(root, log=quiet)                            # net charge +0.5
    with tempfile.TemporaryDirectory() as d:
        root = os.path.join(d, "c")
        _, n = make_system(root, dt=0.002)
        assert not C.check_system(root, log=quiet)                            # wrong production time step
        root = os.path.join(d, "e")
        _, n = make_system(root)
        assert not C.check_system(root, expect_atoms=n + 1, log=quiet)        # atom count mismatch
        assert not C.check_system(root, expect_box=(100, 82, 82), log=quiet)  # box mismatch
        root = os.path.join(d, "g")
        _, n = make_system(root)
        gro = os.path.join(root, "gromacs", "step3_input.gro")
        lines = open(gro).read().splitlines()
        lines[1] = str(n + 5)
        open(gro, "w").write("\n".join(lines) + "\n")
        assert not C.check_system(root, log=quiet)                            # .gro and .psf disagree
        root = os.path.join(d, "f")
        make_system(root)
        os.remove(os.path.join(root, "openmm", "step5_production.inp"))
        open(os.path.join(root, "openmm", "step5_production.inp"), "w").write("nstep = 1\ndt = 0.004\ntemp = 300\n")
        assert not C.check_system(root, log=quiet)                            # settings missing


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print("PASS", t.__name__)
    print(f"{len(tests)} tests passed")
