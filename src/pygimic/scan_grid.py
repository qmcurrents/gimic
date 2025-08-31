#!/usr/bin/env python3

import os, shutil, sys, numpy as np
import copy, time

XYZ2NUM = {"X": 0, "Y": 1, "Z": 2}
OUT_FILES = ("dia_scan.dat", "para_scan.dat", "total_scan.dat")

def build_angle_list(step, max_angle, nsteps):
    """Return a 1‑D array of scan angles (degrees)."""
    if step and nsteps:
        return np.linspace(0.0, step * (nsteps - 1), nsteps)

    if step and max_angle:
        return np.arange(0.0, max_angle + 1e-12, step)

    if nsteps and max_angle:
        return np.linspace(0.0, max_angle, nsteps)

    sys.stderr.write("Need at least two of: step / max_angle / nsteps.\n")
    sys.exit(1)


def collect_groups(xdens_root, add_total, workdir):
    """Return the list of sub‑directories that really contain an XDENS file."""
    if add_total == 'True':
        total_dir = os.path.join(xdens_root, "Total")
        shutil.rmtree(total_dir, ignore_errors=True)
        os.makedirs(total_dir, exist_ok=True)
        os.symlink(os.path.join(workdir, "XDENS"),
                   os.path.join(total_dir, "XDENS"))
        os.symlink(os.path.join(workdir, "MOL"),
                   os.path.join(total_dir, "MOL"))
    groups_num = []
    groups = []
    for name in sorted(os.listdir(xdens_root)):
        gdir = os.path.join(xdens_root, name)
        if os.path.isfile(os.path.join(gdir, "XDENS")):
            if name.split('_')[0].isdigit():
                groups_num.append(name)
            else:
                groups.append(name)
        else:
            print(f"[WARN] skipping \"{gdir}\" no XDENS found.\n")

    groups = sorted(groups_num, key=lambda s: int(s.split('_')[0])) + groups

    if not groups:
        sys.stderr.write(f"No valid sub‑directories in {xdens_root}.\n")
        sys.exit(1)
    return groups

def run_gimic(args, inkeys):
    """Run gimic without writing input and output."""
    if inkeys.getkw('backend')[0] == 'fgimic':
        from fgimic.gimic import GimicDriver
    else:
        from pygimic.pygimic import GimicDriver
    gimic = GimicDriver(args, inkeys)
    return gimic.run_capture()

def grid_scan_dryrun(angles, args, inkeys, orig_grid_sec, init_euler, rot_axis_idx):
    """Make grid grid_scan.xyz containing rotated integration plane."""
    old_flag_dryrun = inkeys.getkw('dryrun')[0]
    inkeys.setkw('dryrun', 'True')

    grid_scan = []
    for j, ang in enumerate(angles):
        grid = copy.deepcopy(orig_grid_sec)
        euler = init_euler.copy()
        euler[rot_axis_idx] += ang
        grid.setkw('rotation', [str(x) for x in euler])

        tmp_inkeys = copy.deepcopy(inkeys)
        tmp_inkeys.add_sect(grid, set=True)

        out_text, _, _ = run_gimic(args, tmp_inkeys)
        if j ==0 and (args.dryrun or old_flag_dryrun == 'True'):
            print(out_text)

        with open('grid.xyz') as f:
            text = f.readlines()
            text[1] = f"Angle= {ang: 8.3f}\n"
            grid_scan.append(text)

    inkeys.setkw('dryrun', old_flag_dryrun)

    with open('grid_scan.xyz', 'w') as f:
        for text in grid_scan:
            for line in text:
                f.write(line)
            f.write("\n")


def trapz_avg(xs, ys):
    """Averaging using simple trapezoidal rule."""
    acc = sum(0.5 * (xs[i] - xs[i - 1]) * (ys[i] + ys[i - 1])
              for i in range(1, len(xs)))
    return acc / abs(xs[-1] - xs[0])


def parse_currents(out_text):
    """Return dia, para, total currents + header lines."""
    lines = out_text.splitlines()
    for idx, ln in enumerate(lines):
        if "Induced current (au)" in ln:
            dia   = float(lines[idx + 1].split()[4])
            para  = float(lines[idx + 2].split()[4])
            total = float(lines[idx + 4].split()[4])
            return dia, para, total, lines[:idx]
    sys.stderr.write("Failed to parse GIMIC output\n")
    sys.exit(1)


def process_group(group, angles, base_grid, init_euler, rot_axis_idx,
                  inkeys_template, args, xdens_root, i):

    n_ang = len(angles)
    contrib = np.zeros((3, n_ang))
    start_group = time.perf_counter()
    progress = max(1, n_ang // 10)

    for j, ang in enumerate(angles):
        grid = copy.deepcopy(base_grid)
        eul  = init_euler.copy()
        eul[rot_axis_idx] += ang
        grid.setkw('rotation', [str(x) for x in eul])

        tmp_kw = copy.deepcopy(inkeys_template)
        tmp_kw.add_sect(grid, set=True)
        tmp_kw.setkw('basis', os.path.join(xdens_root, group, "MOL"))
        tmp_kw.setkw('xdens', os.path.join(xdens_root, group, "XDENS"))

        out_txt, _, _ = run_gimic(args, tmp_kw)
        dia, para, total, top = parse_currents(out_txt)

        # print banner from gimic output for the first gimic calculation
        if i == 0 and j == 0:
            for line in top:
                print(line)

        contrib[:, j] = (dia, para, total)

        done = j + 1
        if done % progress == 0 or done == n_ang:
            elapsed = time.perf_counter() - start_group
            eta = elapsed / done * (n_ang - done) if done < n_ang else 0.0
            perc = done / n_ang
            line = f"[{group:<15}] {done}/{n_ang} ({perc:5.1%})  ETA {eta:6.1f}s"
            # print in terminal
            if sys.stdout.isatty():
                sys.stdout.write("\r" + line)
                sys.stdout.flush()
            # write in file
            else:
                sys.stdout.write(line + "\n")
            sys.stdout.flush()

    elapsed = time.perf_counter() - start_group
    print(f"   ...done in {elapsed:6.1f}s")

    avgs = [trapz_avg(angles, contrib[k]) for k in range(3)]
    print(f"   dia={avgs[0]:.6e} para={avgs[1]:.6e} tot={avgs[2]:.6e}")

    return contrib, avgs


def run(scan_sect, args, inkeys):

    workdir = os.path.dirname(os.path.abspath(args.infile))
    xdens_root = scan_sect.getkw('xdenses_path')[0]

    if not os.path.isdir(xdens_root):
        sys.stderr.write(f"XDENSes directory not found: {xdens_root}\n")
        if scan_sect.getkw('add_total_to_scan')[0] == 'True':
            os.makedirs(os.path.join(xdens_root))
            sys.stderr.write(f"XDENSes directory made in: {xdens_root}\n")
        else:
            sys.stderr.write(f"XDENSes for scanning not found.\n")
            sys.exit(1)

    # 1. collect MO‑group directories
    groups = collect_groups(xdens_root,
                            scan_sect.getkw('add_total_to_scan')[0],
                            workdir)

    # 2. build angle list
    step_kw, max_kw, n_kw = (scan_sect.getkw(k) for k in
                             ('step', 'max_angle', 'nsteps'))
    angles = build_angle_list(float(step_kw[0]) if step_kw else None,
                              float(max_kw[0]) if max_kw else None,
                              int(n_kw[0]) if n_kw else None)

    # 4. prepare the Grid section prototype
    orig_grid = inkeys.sect.pop('Grid')[0] # untouched original
    orig_grid.setkw('rotation_origin', scan_sect.getkw('rotation_origin'))

    init_euler = np.array(orig_grid.getkw('rotation'), dtype=float)
    rot_axis_idx = XYZ2NUM[scan_sect.getkw('rotation_axis')[0].upper()]

    # 5. switch to calc = integral for the scan
    inkeys.setkw('calc', 'integral')

    # Prepare grid_scan.xyz
    grid_scan_dryrun(angles, args, inkeys, orig_grid, init_euler, rot_axis_idx)

    if not (args.dryrun or inkeys.getkw('dryrun')[0] == 'True'):
        # 5.loop for xdenses from angles
        n_ang, n_grp = len(angles), len(groups)
        contributions = np.zeros((3, n_ang, n_grp))
        avgs = np.zeros((n_grp, 3))

        for i, grp in enumerate(groups):
            contrib_i, avgs_i = process_group(grp, angles, orig_grid,
                                              init_euler, rot_axis_idx,
                                              inkeys, args, xdens_root, i)
            contributions[:, :, i] = contrib_i
            avgs[i] = avgs_i

        # 6. Write raw scans to files
        angle_col = angles.reshape(-1, 1)
        header = 'angle ' + ' '.join(groups)

        for idx, fname in enumerate(OUT_FILES):     # dia, para, total
            data = np.column_stack((angle_col, contributions[idx]))
            np.savetxt(fname, data, header=header, fmt="%12.6e")
        print("\nData written:", ", ".join(OUT_FILES))

        # 7. Pint summary table
        print("\nAveraged contributions  (nA/T)")
        print(f"{'Group':<20} {'dia':>12} {'para':>12} {'total':>12}")
        for grp, (d, p, t) in zip(groups, avgs):
            print(f"{grp:<20} {d:12.4e} {p:12.4e} {t:12.4e}")
        print("\nScan finished.")