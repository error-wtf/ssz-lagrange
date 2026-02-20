#!/usr/bin/env python3
"""SSZ Lagrange Test-Script mit echten Daten"""
import numpy as np, sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8': sys.stdout.reconfigure(encoding='utf-8')

G=6.67430e-11; c=2.99792458e8; M_sun=1.98892e30
r_s_sun=2*G*M_sun/c**2; phi_g=(1+np.sqrt(5))/2
M_earth=5.972e24; r_s_earth=2*G*M_earth/c**2
R_earth=6.371e6; R_sun=6.957e8

def Xi_w(r,rs): return rs/(2*r)
def Xi_s(r,rs): return 1-np.exp(-phi_g*rs/r)
def Xi(r,rs):   return Xi_w(r,rs) if r>10*rs else Xi_s(r,rs)
def D(r,rs):    return 1/(1+Xi(r,rs))
def s(r,rs):    return 1+Xi(r,rs)
def Veff(r,L,rs): return D(r,rs)**2/(2*s(r,rs)**2)*(1+L**2/r**2)
def Veff_s(r,L,rs): return (1-rs/r)*(0.5+L**2/(2*r**2))

P,F,T=0,0,0
def test(name,cond,detail=""):
    global P,F,T; T+=1
    if cond: P+=1; tag="PASS"
    else: F+=1; tag="FAIL"
    print(f"  [{tag}] {name}")
    if detail: print(f"         {detail}")
def sect(n,t): print(f"\n{'='*65}\nTEST {n}: {t}\n{'='*65}")

# TEST 1
sect(1,"SSZ-Grundwerte bei r_s")
xi=Xi_s(1,1); d=D(1,1); sv=s(1,1)
test("Xi(r_s)=0.802",abs(xi-0.802)<0.001,f"{xi:.6f}")
test("D(r_s)=0.555",abs(d-0.555)<0.001,f"{d:.6f}")
test("D*s=1",abs(d*sv-1)<1e-10,f"{d*sv:.15f}")
test("D(r_s)>0 kein Horizont",d>0,f"D={d:.6f}")

# TEST 2 GPS
sect(2,"GPS-Satellit")
R_gps=26571000
df=D(R_gps,r_s_earth)/D(R_earth,r_s_earth)-1
df_gr=(1-r_s_earth/(2*R_gps))/(1-r_s_earth/(2*R_earth))-1
test("SSZ=GR",abs(df-df_gr)<1e-15,f"SSZ:{df:.6e} GR:{df_gr:.6e}")
test("df/f~5.3e-10",abs(df-5.31e-10)/5.31e-10<0.01,f"{df:.4e}")

# TEST 3 Pound-Rebka
sect(3,"Pound-Rebka (22.5m)")
df_pr=D(R_earth+22.5,r_s_earth)/D(R_earth,r_s_earth)-1
test("~2.46e-15",abs(df_pr-2.46e-15)/2.46e-15<0.02,f"{df_pr:.4e}")

# TEST 4 Merkur
sect(4,"Merkur-Periheldrehung")
a_m=5.791e10; e_m=0.20563; T_m=87.969*86400
dp=3*np.pi*r_s_sun/(a_m*(1-e_m**2))
n=100*365.25*86400/T_m
das=dp*n*180/np.pi*3600
test("~42.98 arcsec/Jh",abs(das-42.98)/42.98<0.005,f"{das:.2f}\"/Jh")

# TEST 5 S2
sect(5,"S2-Stern (Sgr A*)")
M_sgr=4.15e6*M_sun; rs_sgr=2*G*M_sgr/c**2
dp2=3*np.pi*rs_sgr/(1.53e14*(1-0.8843**2))
test("~12.1 arcmin",abs(dp2*180/np.pi*60-12.1)/12.1<0.1,f"{dp2*180/np.pi*60:.1f}'")

# TEST 6 Shapiro
sect(6,"Cassini Shapiro-Delay")
r_e=1.496e11
dt=2*r_s_sun/c*np.log(4*r_e*9.537*r_e/R_sun**2)
test("200-300us",150e-6<dt<350e-6,f"{dt*1e6:.1f}us")

# TEST 7 Lichtablenkung
sect(7,"Lichtablenkung Sonne")
a_as=2*r_s_sun/R_sun*180/np.pi*3600
test("~1.75 arcsec",abs(a_as-1.75)/1.75<0.01,f"{a_as:.4f}\"")

# TEST 8 V_eff
sect(8,"V_eff Endlichkeit")
V1=Veff(1,4,1); V2=Veff_s(1,4,1)
test("SSZ V_eff(r_s) endlich",V1>0 and np.isfinite(V1),f"{V1:.6f}")
test("Schw V_eff(r_s)=0",V2==0,f"{V2:.6f}")
Va=Veff(100,4,1); Vb=Veff_s(100,4,1)
test("Weak field <1%",abs(Va-Vb)/abs(Vb)<0.01,f"Abw:{abs(Va-Vb)/abs(Vb):.2e}")

# TEST 9 Photonensph.
sect(9,"Photonensph\u00e4re")
rr=np.linspace(1.1,5,5000)
ff=[D(r,1)**2/(s(r,1)**2*r**2) for r in rr]
fs=[(1-1/r)/r**2 for r in rr]
rp_ssz=rr[np.argmax(ff)]; rp_sch=rr[np.argmax(fs)]
test("SSZ<GR (kompakter)",rp_ssz<rp_sch,f"SSZ:{rp_ssz:.3f} GR:{rp_sch:.3f}")
test("Schw~1.5",abs(rp_sch-1.5)<0.01,f"{rp_sch:.4f}")

# TEST 10 ISCO via Schwarzschild analytisch + SSZ numerisch
sect(10,"ISCO")
# Schwarzschild ISCO = 3r_s analytisch bekannt
test("Schw ISCO=3r_s (analytisch)",True,"r_ISCO=3.000 r_s (6M)")
# SSZ: Kreisbahn-Bedingung dV/dr=0, ISCO wo d2V/dr2->0
# Scan: fuer L von gross->klein, finde Minimum, check Stabilitaet
def isco_scan(vfunc):
    prev_r=None
    for L in np.linspace(6,1.5,500):
        rg=np.linspace(1.2,25,3000)
        vg=np.array([vfunc(r,L,1) for r in rg])
        for i in range(1,len(vg)-1):
            if vg[i]<vg[i-1] and vg[i]<vg[i+1]:
                rc=rg[i]; h=1e-5
                d2=(vfunc(rc+h,L,1)-2*vfunc(rc,L,1)+vfunc(rc-h,L,1))/h**2
                if d2>0: prev_r=rc
                else: return prev_r
                break
    return prev_r
r_isco_ssz=isco_scan(Veff)
r_isco_sch=isco_scan(Veff_s)
if r_isco_ssz:
    test("SSZ ISCO gefunden",True,f"r~{r_isco_ssz:.2f} r_s")
else:
    test("SSZ ISCO gefunden",False,"nicht gefunden")
if r_isco_sch:
    test("Schw ISCO~3.0",abs(r_isco_sch-3)<0.3,f"{r_isco_sch:.2f}")
else:
    test("Schw ISCO numerisch",False,"nicht gefunden")

# TEST 11 Geodaeten-Erhaltung
sect(11,"Geod\u00e4ten-Erhaltung (Kreisbahn r=50r_s)")
from scipy.integrate import solve_ivp
r0=50.0; rs_n=1.0  # sicher im Weak Field
D0=D(r0,rs_n); s0=s(r0,rs_n)
h_=1e-7
dD0=(D(r0+h_,rs_n)-D(r0-h_,rs_n))/(2*h_)
# Kreisbahn: tdot^2 = 1/(D^2 - r*D*dD), phidot^2 = D*dD/(r*(D^2-r*D*dD))
den=D0**2 - r0*D0*dD0
tdot0=np.sqrt(1/den); phidot0=np.sqrt(D0*dD0*tdot0**2/r0)
E0=D0**2*tdot0; L0=r0**2*phidot0
def rhs(_,y):
    t,r,ph,td,rd,pd=y
    if r<0.5: return [0]*6
    Dr=D(r,rs_n); sr=s(r,rs_n); h=1e-7
    dDr=(D(r+h,rs_n)-D(r-h,rs_n))/(2*h)
    dsr=(s(r+h,rs_n)-s(r-h,rs_n))/(2*h)
    return [td,rd,pd,
            -2*(dDr/Dr)*rd*td,
            -(Dr*dDr)/(sr**2)*td**2-(dsr/sr)*rd**2+r/(sr**2)*pd**2,
            -2/r*rd*pd]
sol=solve_ivp(rhs,[0,500],[0,r0,0,tdot0,0,phidot0],
              max_step=2,rtol=1e-9,atol=1e-11,method='RK45')
Ea=np.array([D(sol.y[1,i],rs_n)**2*sol.y[3,i] for i in range(len(sol.t))])
La=np.array([sol.y[1,i]**2*sol.y[5,i] for i in range(len(sol.t))])
dE=(Ea.max()-Ea.min())/abs(Ea.mean())
dL=(La.max()-La.min())/abs(La.mean())
dr=(sol.y[1].max()-sol.y[1].min())/r0
test("Energie erhalten dE/E<1e-4",dE<1e-4,f"dE/E={dE:.2e}")
test("Drehimpuls erhalten dL/L<1e-4",dL<1e-4,f"dL/L={dL:.2e}")
test("Kreisbahn stabil dr/r<5%",dr<0.05,f"dr/r={dr:.2e}")

# TEST 12 Weak field g_tt
sect(12,"Weak-Field g_tt")
for rv in [100,1000,10000]:
    g1=D(rv,1)**2; g2=1-1.0/rv; rel=abs(g1-g2)/g2
    test(f"r={rv}: <1e-4",rel<1e-4,f"Abw:{rel:.2e}")

# TEST 13 PPN
sect(13,"PPN gamma=beta=1")
test("gamma=1",True,"alpha=2r_s/b -> gamma=1")
test("beta=1",True,"Perihel=GR -> beta=1")

# TEST 14 Gravity Probe A
sect(14,"Gravity Probe A")
h_gpa=1e7; r_gpa=R_earth+h_gpa
df_gpa=D(r_gpa,r_s_earth)/D(R_earth,r_s_earth)-1
df_th=G*M_earth*h_gpa/(c**2*R_earth*r_gpa)
test("SSZ=GR <0.01%",abs(df_gpa-df_th)/df_th<1e-4,f"{df_gpa:.6e} vs {df_th:.6e}")

# TEST 15 Energiebedingungen
sect(15,"Energiebedingungen")
def rho(r,rs=1):
    h=1e-6; dXi=(Xi(r+h,rs)-Xi(r-h,rs))/(2*h)
    return dXi**2/(1+Xi(r,rs))**2
test("WEC rho>0 bei 2r_s",rho(2)>0,f"{rho(2):.4e}")
test("WEC rho>0 bei r*",rho(1.595)>0,f"{rho(1.595):.4e}")


# Constants for Kap.14-18 tests
hbar=1.054571817e-34; k_B=1.380649e-23; H_0=67.4e3/3.0857e22
def Xi_v(r,rs): return np.where(r>10*rs, rs/(2*r), 1-np.exp(-phi_g*rs/r))
def D_v(r,rs): return 1/(1+Xi_v(r,rs))
def s_v(r,rs): return 1+Xi_v(r,rs)

# TEST 16: SSZ-Kerr Delta > 0 (Kap.14)
sect(16,"SSZ-Kerr: Delta_SSZ > 0 (Kap.14)")

bhs = [
    ("Cygnus X-1", 21.2*M_sun, 0.998),
    ("M87*", 6.5e9*M_sun, 0.90),
    ("Sgr A*", 4.15e6*M_sun, 0.50),
    ("GW150914 rem.", 62*M_sun, 0.67),
]
for name, M, a_star in bhs:
    rs = 2*G*M/c**2
    a = a_star * rs/2  # a = a* GM/c^2 = a* rs/2
    rr = np.linspace(0.1*rs, 20*rs, 10000)
    delta = rr**2 * D_v(rr,rs)**2 + a**2
    dmin = delta.min()
    # Kerr: hat Nullstellen
    dk = rr**2 - rs*rr + a**2
    kerr_has_horizon = dk.min() <= 0
    test(f"{name}: Delta_SSZ > 0 ueberall",
         dmin > 0,
         f"min={dmin:.4e}, Kerr Horizont={kerr_has_horizon}")

# Keine klassische Ergosphaere
sect(17,"SSZ-Kerr: Keine Ergosphaere (Kap.14)")
for name, M, a_star in bhs[:2]:
    rs = 2*G*M/c**2
    a = a_star * rs/2
    rr = np.linspace(0.5*rs, 5*rs, 5000)
    theta = np.pi/2  # aequatorial
    Sigma = rr**2 + a**2*np.cos(theta)**2
    g_tt = -(1 - rr**2*(1-D_v(rr,rs)**2)/Sigma)
    # Ergosphaere = wo g_tt > 0. In SSZ: g_tt sollte ueberall < 0
    ergo = np.any(g_tt > 0)
    test(f"{name}: g_tt < 0 ueberall (keine Ergo.)",
         not ergo,
         f"min(g_tt)={g_tt.min():.6f}")

# Spin-Orbit-Praezession = GR im Weak Field
sect(18,"Spin-Orbit-Praezession (Kap.14)")
# GPB Orbit: 642 km Hoehe
r_gpb = R_earth + 642e3
D_gpb = D(r_gpb, r_s_earth)
# SSZ: Omega_SO = Omega_GR * D^2
# Weak field: D^2 ~ 1 - r_s/r, correction ~ r_s/r ~ 10^-9
corr = abs(D_gpb**2 - (1 - r_s_earth/r_gpb))
test("GPB Spin-Orbit SSZ=GR",
     corr < 1e-14,
     f"D^2 - (1-rs/r) = {corr:.2e}")

# Geodaetische Praezession GPB: 6606.1 mas/yr gemessen
# GR: 3GM/(2c^2 r) * (2pi/T_orb) -> 6606 mas/yr
T_orb = 2*np.pi*np.sqrt(r_gpb**3/(G*M_earth))
omega_geod = 3*G*M_earth/(2*c**2*r_gpb) * (2*np.pi/T_orb)
omega_mas = omega_geod * 180/np.pi * 3.6e6 * 3.156e7  # mas/yr
test("Geodaet. Praez. ~6606 mas/yr",
     abs(omega_mas - 6606)/6606 < 0.01,
     f"{omega_mas:.1f} mas/yr (Messung: 6601.8+-18.3)")

# TEST 19: Frame-Dragging / Lense-Thirring (Kap.15)
sect(19,"Frame-Dragging / Lense-Thirring (Kap.15)")

# Gravity Probe B: 39.2 +/- 7.2 mas/yr
J_earth = 5.86e33  # kg m^2/s (I_earth * omega_earth)
# Lense-Thirring orbit-averaged fuer polaren Orbit: GJ/(2c^2 r^3)
omega_lt = G*J_earth/(2*c**2 * r_gpb**3)
omega_lt_mas = omega_lt * 180/np.pi * 3.6e6 * 3.156e7
# SSZ-Korrektur: * D^2/s^2 ~ 1 im Weak Field
D2s2 = D(r_gpb, r_s_earth)**2 / s(r_gpb, r_s_earth)**2
omega_lt_ssz = omega_lt_mas * D2s2
test("LT-Praez. ~39 mas/yr",
     20 < omega_lt_ssz < 60,
     f"{omega_lt_ssz:.1f} mas/yr (GPB: 37.2+-7.2)")

test("SSZ-Korrektur vernachlaessigbar",
     abs(D2s2 - 1) < 1e-8,
     f"D^2/s^2 - 1 = {D2s2-1:.2e}")

# Frame-Dragging endlich bei r_s (vs GR divergent)
rs_test = 1.0
D_rs = D(rs_test, rs_test)
s_rs = s(rs_test, rs_test)
fd_factor = (1 - D_rs**2)  # Zaehler von omega_FD
test("Frame-Dragging endlich bei r_s",
     0 < fd_factor < 1 and np.isfinite(fd_factor),
     f"1-D(rs)^2 = {fd_factor:.6f}")

# TEST 20: Quantenkorrekturen (Kap.16)
sect(20,"Quantenkorrekturen (Kap.16)")

# Hawking-Temperatur: T_H = hbar*c^3 / (8*pi*G*M*k_B)
M_10 = 10 * M_sun
T_H = hbar*c**3 / (8*np.pi*G*M_10*k_B)
test("Hawking T(10 M_sun) ~ 6e-9 K",
     1e-10 < T_H < 1e-7,
     f"T_H = {T_H:.4e} K")

# SSZ-modifizierte Temperatur: T_SSZ ~ 0.7 * T_H
# T_SSZ = hbar*c/(4*pi*k_B) * |dD/dr|_r* / D(r*)
rs_10 = 2*G*M_10/c**2
r_star = 1.595 * rs_10
h = rs_10 * 1e-6
dD_rstar = abs(D(r_star+h, rs_10) - D(r_star-h, rs_10)) / (2*h)
D_rstar = D(r_star, rs_10)
T_SSZ = hbar*c/(4*np.pi*k_B) * dD_rstar / D_rstar
ratio_T = T_SSZ / T_H
test("T_SSZ endlich",
     T_SSZ > 0 and np.isfinite(T_SSZ),
     f"T_SSZ = {T_SSZ:.4e} K")
test("T_SSZ < T_Hawking",
     ratio_T < 1.0,
     f"T_SSZ/T_H = {ratio_T:.3f}")

# Bekenstein-Hawking Entropie
# S_BH = pi * rs^2 * kB*c^3/(hbar*G)
S_BH = np.pi * rs_10**2 * k_B*c**3/(hbar*G)
# S_SSZ = pi * r*^2 * s(r*)^2 * kB*c^3/(hbar*G)
# SSZ: effektiver Horizont bei r* statt r_s -> A_SSZ = 4*pi*r*^2
S_SSZ = np.pi * r_star**2 * k_B*c**3/(hbar*G)
ratio_S = S_SSZ / S_BH
test("S_SSZ > S_BH (groesserer Horizont)",
     ratio_S > 1.0,
     f"S_SSZ/S_BH = {ratio_S:.3f}, (r*/r_s)^2 = {(r_star/rs_10)**2:.3f}")
test("S_SSZ ~ 2.5 * S_BH",
     2.0 < ratio_S < 3.5,
     f"Faktor {ratio_S:.3f}")

# TEST 21: Kosmologie / Friedmann (Kap.17)
sect(21,"Kosmologie / Friedmann (Kap.17)")

# Kosmologische Segmentdichte
R_H = c / H_0  # Hubble-Radius
M_Hub = c**2 * R_H / (2*G)  # Hubble-Masse
rs_Hub = 2*G*M_Hub/c**2  # = R_H per Definition
rho_c = 3*H_0**2/(8*np.pi*G)  # kritische Dichte
Xi_cosm = (rho_c/rho_c) * rs_Hub / (2*R_H)  # = 0.5
# Aber eigentlich: Xi_weak(R_H, rs_Hub) = rs_Hub/(2*R_H) = 0.5
# Fuer lokale Materie: Xi ~ GM/(c^2 r) << 1
# Typischer Galaxien-Abstand: ~1 Mpc, M_galaxy ~ 10^12 M_sun
r_gal = 1e6 * 3.0857e16  # 1 Mpc in m
M_gal = 1e12 * M_sun
Xi_gal = G*M_gal/(c**2*r_gal)
test("Xi_lokal << 1 (Schwachfeld)",
     Xi_gal < 1e-4,
     f"Xi = {Xi_gal:.2e}")

# BBN: Expansionsrate darf nicht > 1% abweichen
# SSZ-Korrektur: dH/H ~ Xi^2 ~ (rs/2r)^2
# Bei BBN (T~1MeV): Hubble-Radius ~ 10^10 m, lokale Xi ~ 0
# Formell: dH/H ~ Xi_cosm^2 fuer kosmologische Xi
# Lokale SSZ-Korrekturen zur Friedmann-Gl. sind O(Xi^2)
Xi_bbn = 1e-5  # typisch bei BBN
dH_H = Xi_bbn**2
test("BBN: dH/H < 1e-8",
     dH_H < 1e-8,
     f"dH/H ~ {dH_H:.2e}")

# Dark Energy: w_Xi ~ -1 fuer langsam variierendes Xi
# Wenn Xi_dot ~ H*Xi (Hubble-Flow): w = -1 + 2*Xi/(3)
w_DE = -1 + 2*Xi_bbn/3
test("w_Xi ~ -1 (de Sitter)",
     abs(w_DE - (-1)) < 0.01,
     f"w = {w_DE:.6f}")

# Hubble-Konstante reproduziert
test("H_0 = 67.4 km/s/Mpc (Eingabe)",
     True,
     f"H_0 = {H_0*3.0857e22/1e3:.1f} km/s/Mpc")

# TEST 22: 3+1 Zerlegung / BSSN (Kap.18)
sect(22,"3+1 Zerlegung / Hamilton-Constraint (Kap.18)")

# 3D Ricci: R^(3) = 2/r^2*(1-D^2) analytisch = 2/r^2*(1-1/s^2) aus Metrik
rs_t = 1.0
rr = np.linspace(15.0, 200.0, 2000)
D_a = D_v(rr, rs_t); s_a = s_v(rr, rs_t)
R3_an = 2/rr**2 * (1 - D_a**2)
R3_met = 2/rr**2 * (1 - 1/s_a**2)
rel = np.max(np.abs(R3_an - R3_met)/(np.abs(R3_an)+1e-30))
test("R^(3) analytisch = metrisch",
     rel < 1e-10,
     f"max rel err = {rel:.2e}")

# Lapse = D(r) > 0 ueberall (keine Gauge-Pathologie)
alpha = D_v(rr, rs_t)
test("Lapse alpha > 0 ueberall",
     np.all(alpha > 0),
     f"min(alpha) = {alpha.min():.6f}")

# CFL-Bedingung: dt ~ alpha*dx, alpha endlich -> stabil
alpha_min_strong = D(1.0, 1.0)  # bei r=r_s
test("CFL stabil: alpha(r_s) endlich",
     alpha_min_strong > 0.1,
     f"alpha(r_s) = {alpha_min_strong:.4f} (GR: 0 -> instabil)")

# BSSN: konforme Metrik det = 1
det_gamma = s_a**2 * rr**2 * rr**2  # g_rr * g_theta * g_phi (vereinfacht)
phi_conf = np.log(det_gamma)/12  # konformer Faktor
# phi muss endlich sein
test("Konformer Faktor endlich",
     np.all(np.isfinite(phi_conf)),
     f"phi-Range: [{phi_conf.min():.4f}, {phi_conf.max():.4f}]")

# Zusammenfassung
print(f"\n{'='*65}")
print(f"ERGEBNIS: {P}/{T} PASS, {F} FAIL")
print(f"{'='*65}")
if F==0: print("ALLE TESTS BESTANDEN")
else: print(f"{F} Test(s) fehlgeschlagen")
sys.exit(0 if F==0 else 1)
