#!/usr/bin/env python3
"""Run SSZ cosmology scaffold and produce results."""
import ssz_cosmo as cosmo
import json, csv, numpy as np

def main():
    zs = [0, 0.5, 1, 2, 5, 10]
    results = {}
    for mode in ['divide', 'multiply']:
        rows = []
        for z in zs:
            D_gr, f_gr = cosmo.growth_factor_gr(z)
            D_ssz, f_ssz = cosmo.growth_factor_ssz(z, mode=mode)
            Hr = cosmo.H_ratio(z, mode=mode)
            rows.append({'z': z, 'D_GR': D_gr, 'D_SSZ': D_ssz,
                         'f_GR': f_gr, 'f_SSZ': f_ssz,
                         'Hratio_SSZ/GR': Hr})
        results[mode] = rows
        fn = f'results/ssz_cosmo_results_{mode}.csv'
        with open(fn, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=rows[0].keys())
            w.writeheader(); w.writerows(rows)
        print(f'Wrote {fn}')

    # Summary
    summary = {}
    for mode in ['divide', 'multiply']:
        cmb = cosmo.cmb_acoustic_proxy(mode=mode)
        bbn = cosmo.bbn_proxy(mode=mode)
        summary[mode] = {'cmb': cmb, 'bbn': bbn}
    with open('results/ssz_cosmo_summary_all.json', 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print('Wrote results/ssz_cosmo_summary_all.json')

    # Plots
    try:
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        for mode in ['divide', 'multiply']:
            zp = np.linspace(0, 10, 200)
            Hr = [cosmo.H_ratio(z, mode=mode) for z in zp]
            ax1.plot(zp, Hr, label=f'{mode}')
        ax1.set_xlabel('z'); ax1.set_ylabel('H_SSZ/H_GR')
        ax1.set_title('Hubble ratio'); ax1.legend()
        ax1.axhline(1, ls='--', c='gray')

        for mode in ['divide', 'multiply']:
            zp = np.linspace(0, 10, 200)
            Dg = [cosmo.growth_factor_ssz(z, mode=mode)[0] for z in zp]
            Dgr = [cosmo.growth_factor_gr(z)[0] for z in zp]
            ax2.plot(zp, Dg, label=f'SSZ {mode}')
        ax2.plot(zp, Dgr, 'k--', label='GR')
        ax2.set_xlabel('z'); ax2.set_ylabel('D(z)')
        ax2.set_title('Growth factor'); ax2.legend()

        plt.tight_layout()
        plt.savefig('results/growth_and_Hratio.png', dpi=150)
        print('Wrote results/growth_and_Hratio.png')
    except ImportError:
        print('matplotlib not available, skipping plots')

if __name__ == '__main__':
    main()
