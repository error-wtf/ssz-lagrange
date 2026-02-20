
import numpy as np
import pandas as pd
from ssz_cosmo import CosmoParams, acoustic_scale_theta, bbn_H_ratio, growth_factor_LCDM, growth_rate_f, E_GR, E_SSZ

def run_one(mode: str):
    p = CosmoParams(coupling_mode=mode)  # defaults + mode
    z_star = 1090.0

    gr = acoustic_scale_theta(z_star, p, use_ssz_H=False)
    ssz = acoustic_scale_theta(z_star, p, use_ssz_H=True)
    bbn_ratio = bbn_H_ratio(p)

    zs = np.array([0, 0.5, 1, 2, 5, 10], dtype=float)
    D_gr = growth_factor_LCDM(zs, p, use_ssz_H=False)
    D_ssz = growth_factor_LCDM(zs, p, use_ssz_H=True)
    f_gr = growth_rate_f(zs, p, use_ssz_H=False)
    f_ssz = growth_rate_f(zs, p, use_ssz_H=True)

    df = pd.DataFrame({
        "z": zs,
        "D_GR": D_gr,
        "D_SSZ": D_ssz,
        "f_GR": f_gr,
        "f_SSZ": f_ssz,
        "Hratio_SSZ/GR": E_SSZ(zs, p)/E_GR(zs, p),
    })

    summary = {
        "coupling_mode": mode,
        "params": p.__dict__,
        "CMB_theta_GR_rad": gr["theta_rad"],
        "CMB_theta_SSZ_rad": ssz["theta_rad"],
        "CMB_theta_ratio": ssz["theta_rad"]/gr["theta_rad"],
        "CMB_rs_ratio": ssz["r_s_m"]/gr["r_s_m"],
        "CMB_DA_ratio": ssz["D_A_m"]/gr["D_A_m"],
        "BBN_H_ratio_at_z=4e9": bbn_ratio,
    }
    return df, summary

def main():
    outs = {}
    for mode in ["divide", "multiply"]:
        df, summary = run_one(mode)
        df.to_csv(f"ssz_cosmo_results_{mode}.csv", index=False)
        with open(f"ssz_cosmo_summary_{mode}.json","w",encoding="utf-8") as f:
            import json; json.dump(summary, f, indent=2)
        outs[mode] = summary

    # write a combined summary
    with open("ssz_cosmo_summary_all.json","w",encoding="utf-8") as f:
        import json; json.dump(outs, f, indent=2)

    print("=== SSZ cosmology scaffold (both coupling modes) ===")
    print(outs)

if __name__ == "__main__":
    main()
