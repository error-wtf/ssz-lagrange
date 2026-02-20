
"""
SSZ Cosmology prototype (agent-mode scaffold)

Goal:
- Provide a *testable* pipeline for:
  1) CMB acoustic-scale signatures (theta_*, r_s, D_A)
  2) BBN expansion-rate constraint (H at ~1 MeV epoch proxy via very high z)
  3) Linear structure growth (growth factor D1(a), growth rate f)
  4) "No free parameters" check: SSZ adds no new fit parameters beyond standard cosmology inputs.
     Any SSZ mapping must be fixed by SSZ theory; here we provide a *minimal* mapping:
        Xi_cosmo(z) ~= |Phi(z)|/c^2,
     motivated by SSZ weak-field identity Xi ~ GM/(c^2 r) ~ |Phi|/c^2.

This module intentionally avoids curve-fitting. It is a forward model:
you pick parameters and compute residuals against data.

Default cosmological parameters are provided as common reference values (editable by user).
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np

C = 299792458.0  # m/s
G = 6.67430e-11  # SI

# --- Radiation density for photons from Tcmb; neutrinos approximated with Neff if needed.
# We compute Omega_gamma h^2 from Tcmb, using standard formula.
# This is "standard cosmology", not SSZ-specific.

@dataclass(frozen=True)
class CosmoParams:
    h: float = 0.674
    Omega_m: float = 0.315
    Omega_b: float = 0.049
    Omega_L: float = 1.0 - 0.315  # flat by default
    Tcmb: float = 2.7255
    Neff: float = 3.046
    # Scalar potential RMS amplitude proxy (dimensionless Phi/c^2)
    # Using ~1e-5 motivated by CMB anisotropy scale. This is NOT a fit; change if you want to test.
    Phi0: float = 1.0e-5
    # Switch: if set, overrides Xi(z) with constant value (useful falsification test)
    Xi_const: float | None = None

    # How SSZ time dilation couples to background expansion:
    # - 'divide'   : H_SSZ = H_GR / D(z)  (H defined wrt proper time; da/dtau = (da/dt)/D)
    # - 'multiply' : H_SSZ = H_GR * D(z)  (H defined wrt coordinate time; useful alternative convention)
    coupling_mode: str = 'divide'

    @property
    def H0(self) -> float:
        return 100.0 * self.h * 1000.0 / (3.085677581491367e22)  # 1/s

    @property
    def Omega_r(self) -> float:
        # Omega_gamma h^2 = 2.469e-5 (Tcmb/2.7255)^4
        Omega_gamma_h2 = 2.469e-5 * (self.Tcmb / 2.7255) ** 4
        Omega_gamma = Omega_gamma_h2 / (self.h ** 2)
        # neutrino factor (standard): rho_nu/rho_gamma = (7/8)*(4/11)^(4/3)*Neff
        f_nu = (7.0/8.0) * (4.0/11.0) ** (4.0/3.0) * self.Neff
        return Omega_gamma * (1.0 + f_nu)

def E_GR(z: np.ndarray, p: CosmoParams) -> np.ndarray:
    """Dimensionless H(z)/H0 in flat LCDM (GR background)."""
    zp1 = 1.0 + z
    return np.sqrt(p.Omega_r * zp1**4 + p.Omega_m * zp1**3 + p.Omega_L)

def Xi_cosmo(z: np.ndarray, p: CosmoParams) -> np.ndarray:
    """
    Minimal SSZ cosmology mapping:
      Xi ~ |Phi|/c^2, where Phi is the (dimensionless) Newtonian potential amplitude.
    We model Phi(z) with a very simple heuristic:
      - During matter domination, Phi ~ const.
      - During Lambda domination, Phi decays mildly; we approximate with growth suppression g(z)/a.

    We implement a standard growth-based estimate without introducing new SSZ parameters:
      Phi(z) = Phi0 * D1(z)/a  /  (D1(0)/1)
    i.e., normalized to Phi0 at z=0.
    """
    if p.Xi_const is not None:
        return np.full_like(z, float(p.Xi_const), dtype=float)

    # Compute growth factor D1(z) using GR background as baseline for Phi evolution.
    # This does not fit anything; it's internal consistency for the Phi(z) ansatz.
    # For z up to ~1e4 this is fine; for extreme z we clamp to matter-era behavior.
    z = np.asarray(z, dtype=float)
    if np.max(z) == 0.0:
        return np.ones_like(z, dtype=float)
    # If very high z, Phi ~ const => Xi ~ Phi0 (since D1 ~ a)
    out = np.empty_like(z, dtype=float)

    # Moderate range: compute via ODE solution with GR H(z)
    mask = z <= 2e3
    if np.any(mask):
        zz = z[mask]
        D = growth_factor_LCDM(zz, p, use_ssz_H=False)  # GR growth
        a = 1.0 / (1.0 + zz)
        D0 = growth_factor_LCDM(np.array([0.0]), p, use_ssz_H=False)[0]
        Phi = p.Phi0 * (D / a) / (D0 / 1.0)
        out[mask] = np.abs(Phi)

    # Very high z: matter+radiation era, Phi approximately constant on large scales.
    if np.any(~mask):
        out[~mask] = np.abs(p.Phi0)

    return out

def D_ssz(z: np.ndarray, p: CosmoParams) -> np.ndarray:
    """SSZ time-dilation factor in cosmology: D=1/(1+Xi(z))."""
    xi = Xi_cosmo(z, p)
    return 1.0 / (1.0 + xi)

def E_SSZ(z: np.ndarray, p: CosmoParams) -> np.ndarray:
    """Dimensionless H(z)/H0 under a minimal SSZ coupling.

    IMPORTANT: In cosmology one must be explicit about the time variable:
    - If H is defined with respect to *proper time* of comoving observers, and SSZ says d\u03c4 = D(z) dt,
      then da/d\u03c4 = (da/dt)/D \u21d2 H_SSZ = H_GR / D.
    - If H is instead defined with respect to a chosen *coordinate time* t and you want the proper-time
      effect to appear as a slow-down in H, you may use H_SSZ = H_GR * D.

    This code supports both via p.coupling_mode \u2208 {'divide','multiply'}.
    """
    if p.coupling_mode == 'divide':
        return E_GR(z, p) / D_ssz(z, p)
    if p.coupling_mode == 'multiply':
        return E_GR(z, p) * D_ssz(z, p)
    raise ValueError(f"Unknown coupling_mode={p.coupling_mode!r}")

def comoving_distance(z: float, p: CosmoParams, use_ssz_H: bool = True, n: int = 4000) -> float:
    """Comoving line-of-sight distance chi(z) = c \u222b dz/H(z)."""
    zz = np.linspace(0.0, z, n)
    if use_ssz_H:
        Ez = E_SSZ(zz, p)
    else:
        Ez = E_GR(zz, p)
    integrand = 1.0 / Ez
    return C / p.H0 * np.trapz(integrand, zz)

def angular_diameter_distance(z: float, p: CosmoParams, use_ssz_H: bool = True) -> float:
    """Flat universe: D_A = chi/(1+z)."""
    return comoving_distance(z, p, use_ssz_H=use_ssz_H) / (1.0 + z)

def sound_horizon(z_star: float, p: CosmoParams, use_ssz_H: bool = True, n: int = 6000) -> float:
    """
    Comoving sound horizon:
      r_s(z*) = \u222b_{z*}^{\u221e} c_s(z)/H(z) dz
    We truncate at z_max where radiation dominates and contribution is negligible.
    """
    z_max = 1.0e7
    zz = np.logspace(np.log10(z_star), np.log10(z_max), n)
    # baryon-to-photon ratio R(z) = 3 rho_b / (4 rho_gamma)
    Omega_gamma_h2 = 2.469e-5 * (p.Tcmb / 2.7255) ** 4
    Omega_gamma = Omega_gamma_h2 / (p.h ** 2)
    R = (3.0 * p.Omega_b) / (4.0 * Omega_gamma) * 1.0 / (1.0 + zz)
    c_s = C / np.sqrt(3.0 * (1.0 + R))
    if use_ssz_H:
        Ez = E_SSZ(zz, p)
    else:
        Ez = E_GR(zz, p)
    Hz = p.H0 * Ez
    integrand = c_s / Hz
    return np.trapz(integrand, zz)

def acoustic_scale_theta(z_star: float, p: CosmoParams, use_ssz_H: bool = True) -> dict:
    """Return theta_* = r_s(z*)/D_A(z*), plus components."""
    rs = sound_horizon(z_star, p, use_ssz_H=use_ssz_H)
    da = angular_diameter_distance(z_star, p, use_ssz_H=use_ssz_H)
    return {"z_star": z_star, "r_s_m": rs, "D_A_m": da, "theta_rad": rs/da}

def growth_factor_LCDM(z: np.ndarray, p: CosmoParams, use_ssz_H: bool = True) -> np.ndarray:
    """
    Linear growth factor D1(z) (unnormalized) by solving ODE in ln a:
      d2D/dx2 + (2 + d ln H / dx) dD/dx - (3/2) Omega_m(a) D = 0
    where x = ln a.

    Returns D normalized such that D(0)=1.
    """
    z = np.asarray(z, dtype=float)
    if np.max(z) == 0.0:
        return np.ones_like(z, dtype=float)
    # Solve on a grid in a, then interpolate.
    a_min = 1.0 / (1.0 + np.max(z))
    a_grid = np.logspace(np.log10(a_min), 0.0, 3000)
    x = np.log(a_grid)

    def E_of_a(a):
        zz = 1.0/a - 1.0
        if use_ssz_H:
            return E_SSZ(zz, p)
        else:
            return E_GR(zz, p)

    E = E_of_a(a_grid)
    # d ln H / d ln a = d ln E / d ln a
    dlnE_dx = np.gradient(np.log(E), x)

    # Omega_m(a) = Omega_m a^-3 / E(a)^2
    Om_a = p.Omega_m * a_grid**(-3) / (E**2)

    # ODE as first-order system: y0=D, y1=dD/dx
    y0 = np.zeros_like(a_grid)
    y1 = np.zeros_like(a_grid)

    # Initial conditions in deep matter era: D ~ a, so dD/dx = D
    y0[0] = a_grid[0]
    y1[0] = y0[0]

    dx = np.diff(x)
    for i in range(len(dx)):
        # simple RK4
        def f(state, idx):
            D, dD = state
            A = 2.0 + dlnE_dx[idx]
            B = 1.5 * Om_a[idx]
            # dD/dx = dD
            # ddD/dx = -A dD + B D
            return np.array([dD, -A*dD + B*D])

        s = np.array([y0[i], y1[i]])
        h = dx[i]
        k1 = f(s, i)
        k2 = f(s + 0.5*h*k1, i)
        k3 = f(s + 0.5*h*k2, i)
        k4 = f(s + h*k3, i)
        s_next = s + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        y0[i+1], y1[i+1] = s_next

    # Normalize so D(a=1)=1
    D_norm = y0 / y0[-1]

    # Interpolate to requested z
    a_req = 1.0 / (1.0 + z)
    return np.interp(a_req, a_grid, D_norm)

def growth_rate_f(z: np.ndarray, p: CosmoParams, use_ssz_H: bool = True) -> np.ndarray:
    """f = d ln D / d ln a evaluated numerically."""
    z = np.asarray(z, dtype=float)
    # small step in ln a
    a = 1.0/(1.0+z)
    eps = 1e-4
    z1 = 1.0/(a*np.exp(eps)) - 1.0
    z2 = 1.0/(a*np.exp(-eps)) - 1.0
    D1 = growth_factor_LCDM(z1, p, use_ssz_H=use_ssz_H)
    D2 = growth_factor_LCDM(z2, p, use_ssz_H=use_ssz_H)
    return (np.log(D1) - np.log(D2)) / (2*eps)

def bbn_H_ratio(p: CosmoParams, z_bbn: float = 4.0e9) -> float:
    """
    Proxy for BBN epoch expansion rate ratio H_SSZ/H_GR at a very high redshift.
    Real mapping from T to z depends on entropy, g*, etc. This is *constraint scaffolding*:
    if SSZ modifies H(z) significantly at early times, it will impact BBN yields.
    """
    z = np.array([z_bbn], dtype=float)
    return float(E_SSZ(z, p)[0] / E_GR(z, p)[0])
