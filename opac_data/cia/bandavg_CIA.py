"""Band-averaged CIA tables (form 3) on a model's wavelength bands.

Each band's value is the mean of the HITRAN CIA over the band in wavenumber,
  k_b(T) = 1/(wn_hi - wn_lo) * int_{wn_lo}^{wn_hi} k(wn, T) dwn,
integrated with the trapezoid rule on the HITRAN 1 cm-1 grid (k linear in wn at the band
edges); k = 0 outside the HITRAN wavenumber range, so a band only partly covered by the data
gets the mean over the whole band. (The older reformat_CIA.py sampled k at the band-centre
wavenumber, which differs from the band mean by up to a factor of ~10 for the 11-band grid.)

Output format (read by ck_opacity_mod, CIA form 3, in Exo-FMS_column_ck and Exo-FMS_2):
  nT nb
  T [K] (nT values)
  band-centre wavenumbers [cm-1], ascending (nb values)
  (blank line)
  nT rows of nb values: k [cm5 molecule-2], ascending wavenumber; zero written as 1e-99

Usage: python bandavg_CIA.py ../wavelengths_GCM_11.txt 11
  (run in this directory, with the HITRAN *_2011.cia files present)
"""
import sys
import numpy as np

SPECIES = ['H2-H2', 'H2-He', 'H2-H', 'He-H']
K_FLOOR = 1.0e-99


def read_hitran(fname):
    """Read a HITRAN CIA file: list of (T, wn[nrec], k[nrec]) blocks."""
    blocks = []
    with open(fname) as f:
        while True:
            h = f.readline().split()
            if not h:
                break
            nrec = int(h[3])
            T = float(h[4])
            d = np.array([f.readline().split()[:2] for _ in range(nrec)], dtype=float)
            blocks.append((T, d[:, 0], d[:, 1]))
    return blocks


def band_mean(wn, k, lo, hi):
    """Mean of k over [lo, hi] in wavenumber, k = 0 outside [wn[0], wn[-1]]."""
    a, b = max(lo, wn[0]), min(hi, wn[-1])
    if b <= a:
        return 0.0
    inner = (wn > a) & (wn < b)
    x = np.concatenate(([a], wn[inner], [b]))
    y = np.concatenate(([np.interp(a, wn, k)], k[inner], [np.interp(b, wn, k)]))
    return np.trapezoid(y, x) / (hi - lo)


def main():
    wl_file, tag = sys.argv[1], sys.argv[2]
    wl = np.loadtxt(wl_file, skiprows=1)            # band edges [um], ascending
    wn_e = np.sort(1.0e4 / wl)                        # band edges [cm-1], ascending
    nb = len(wn_e) - 1
    wn_c = 0.5 * (wn_e[:-1] + wn_e[1:])

    for sp in SPECIES:
        blocks = read_hitran(sp + '_2011.cia')
        T = np.array([blk[0] for blk in blocks])
        kb = np.zeros((len(T), nb))
        for i, (_, wn, k) in enumerate(blocks):
            for b in range(nb):
                kb[i, b] = band_mean(wn, k, wn_e[b], wn_e[b + 1])
        kb = np.maximum(kb, K_FLOOR)

        fname = f'{sp}_reform_{tag}.txt'
        with open(fname, 'w') as f:
            f.write(f'{len(T)} {nb} \n')
            f.write(' '.join(str(v) for v in T) + '\n')
            f.write(' '.join(str(v) for v in wn_c) + '\n')
            f.write('\n')
            for i in range(len(T)):
                f.write(' '.join(str(v) for v in kb[i]) + '\n')
        print(f'{fname}: nT = {len(T)}, nb = {nb}')


if __name__ == '__main__':
    main()
