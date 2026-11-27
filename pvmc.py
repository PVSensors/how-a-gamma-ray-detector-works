"""pvmc – a small Monte-Carlo of gamma rays in a cylindrical detector crystal.

Python version of the simulation that runs in the browser on pvsensors.com/detectors
(src/lib/transport.ts), generalised to any detector material and source.

Physics
  * cross-sections: NIST XCOM (photoelectric, incoherent = Compton, pair production),
    compounds by the mixture rule, from data/xcom_elements.json
  * Compton angle and energy sampled from the Klein–Nishina formula (Kahn's method)
  * photoelectric absorption followed by K-shell fluorescence (escape peaks)
  * pair production with two 511 keV annihilation photons (single / double escape)
  * thin aluminium entrance window; optional lead shield that sends back
    backscattered photons, Pb K X-rays and 511 keV annihilation photons
  * Gaussian energy resolution  FWHM(E) = R662 · sqrt(662 keV / E) · E   (with a floor)

Not modelled (keep in mind for semiconductors): charge trapping / low-energy tails,
electron escape and bremsstrahlung, true coincidence summing of cascade gammas,
pile-up, the detailed shape of the shield. It is a teaching model, not Geant4.

© PV Sensors – MIT licence. Data: NIST XCOM, doi:10.18434/T48G6X.
"""
from __future__ import annotations

import json
import math
import random
import re
from dataclasses import dataclass, field
from pathlib import Path

ME_C2 = 511.0  # keV

# --------------------------------------------------------------------------- data
AW = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "Na": 22.990, "Si": 28.085, "S": 32.06,
      "Zn": 65.38, "Ge": 72.630, "Se": 78.971, "Br": 79.904, "Y": 88.906, "Cd": 112.41, "Te": 127.60,
      "I": 126.90, "Cs": 132.91, "La": 138.91, "Gd": 157.25, "Lu": 174.97, "Pb": 207.2, "Bi": 208.98}

# K-shell data used for fluorescence escape: (K-edge keV, mean Kα energy keV, fluorescence yield ωK)
K_SHELL = {"Zn": (9.659, 8.62, 0.47), "Ge": (11.103, 9.87, 0.54), "Se": (12.658, 11.20, 0.59),
           "Br": (13.474, 11.90, 0.62), "Y": (17.038, 14.93, 0.71), "Cd": (26.711, 23.11, 0.84),
           "Te": (31.814, 27.37, 0.875), "I": (33.169, 28.54, 0.884), "Cs": (35.985, 30.85, 0.895),
           "La": (38.925, 33.30, 0.905), "Gd": (50.239, 42.85, 0.935), "Lu": (63.314, 53.92, 0.95),
           "Pb": (88.005, 74.82, 0.963), "Bi": (90.526, 76.86, 0.965)}
K_FRACTION = 0.86   # share of photo-absorption on an element that happens in its K shell (above the edge)


def load_xcom(path="data/xcom_elements.json"):
    return json.loads(Path(path).read_text())


def mass_fractions(formula: str) -> dict:
    parts = re.findall(r"([A-Z][a-z]?)([0-9.]*)", formula)
    m = [(el, (float(n) if n else 1.0) * AW[el]) for el, n in parts]
    tot = sum(x for _, x in m)
    w: dict = {}
    for el, x in m:
        w[el] = w.get(el, 0.0) + x / tot
    return w


# --------------------------------------------------------------------------- materials
@dataclass
class Material:
    """A detector material. Cross-sections are built from the XCOM element tables."""
    name: str
    formula: str
    density: float                 # g/cm³
    fwhm_662: float = 7.0          # energy resolution at 662 keV, % FWHM
    kind: str = "scintillator"     # 'scintillator' or 'semiconductor' (only sets the resolution floor)
    xcom: dict = field(default=None, repr=False)

    def __post_init__(self):
        if self.xcom is None:
            raise ValueError("pass xcom=load_xcom(...)")
        X = self.xcom
        self.E = X["energy_keV"]
        self.w = mass_fractions(self.formula)
        missing = [el for el in self.w if el not in X["elements"]]
        if missing:
            raise ValueError(f"no XCOM data for {missing}; download them with fetch_xcom() in the notebook")
        n = len(self.E)
        col = lambda el, c: X["elements"][el][c]
        self._photo_el = {el: [self.w[el] * col(el, "photoelectric")[i] for i in range(n)] for el in self.w}
        self.photo = [sum(self._photo_el[el][i] for el in self.w) for i in range(n)]
        self.incoh = [sum(self.w[el] * col(el, "incoherent")[i] for el in self.w) for i in range(n)]
        self.pair = [sum(self.w[el] * (col(el, "pair_nuclear")[i] + col(el, "pair_electron")[i]) for el in self.w)
                     for i in range(n)]
        self._logE = [math.log(e) for e in self.E]
        self.k_lines = [(el, *K_SHELL[el]) for el in self.w if el in K_SHELL]

    # log-log interpolation (linear where a value is zero)
    def _at(self, arr, E):
        Es = self.E
        if E <= Es[0]:
            return arr[0]
        lo, hi = 0, len(Es) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if Es[mid] <= E: lo = mid
            else: hi = mid
        a, b = arr[lo], arr[hi]
        if a <= 0 or b <= 0:
            return a + (E - Es[lo]) / (Es[hi] - Es[lo]) * (b - a)
        f = (math.log(E) - self._logE[lo]) / (self._logE[hi] - self._logE[lo])
        return math.exp(math.log(a) + f * (math.log(b) - math.log(a)))

    def mu(self, E):
        """linear attenuation coefficients (1/cm): photo, Compton, pair"""
        r = self.density
        return self._at(self.photo, E) * r, self._at(self.incoh, E) * r, self._at(self.pair, E) * r

    def fluorescence(self, E, rnd):
        """after a photo-absorption at energy E: return the K X-ray energy emitted, or 0"""
        if not self.k_lines:
            return 0.0
        tot = self._at(self.photo, E)
        u = rnd() * tot
        acc = 0.0
        for el, edge, ka, omega in self.k_lines:
            acc += self._at(self._photo_el[el], E)
            if u <= acc:
                if E > edge and rnd() < K_FRACTION * omega:
                    return ka
                return 0.0
        return 0.0

    def fwhm(self, E):
        """FWHM in keV at energy E (same model as pvsensors.com)"""
        rel = self.fwhm_662 / 100 * math.sqrt(662 / max(E, 20))
        floor = 0.8 if self.kind == "semiconductor" else 2.0
        return max(rel * E, floor)


# Catalogue of the crystals on pvsensors.com (density g/cm³, typical resolution at 662 keV).
# Change fwhm_662 to the value of YOUR crystal.
CATALOGUE = {
    "NaI(Tl)":   dict(formula="NaI", density=3.67, fwhm_662=7.0, kind="scintillator"),
    "CsI(Tl)":   dict(formula="CsI", density=4.51, fwhm_662=6.5, kind="scintillator"),
    "LaBr3(Ce)": dict(formula="LaBr3", density=5.08, fwhm_662=2.9, kind="scintillator"),
    "BGO":       dict(formula="Bi4Ge3O12", density=7.13, fwhm_662=10.5, kind="scintillator"),
    "LYSO(Ce)":  dict(formula="Lu1.8Y0.2SiO5", density=7.1, fwhm_662=9.0, kind="scintillator"),
    "GOS":       dict(formula="Gd2O2S", density=7.34, fwhm_662=15.0, kind="scintillator"),
    "PVT plastic": dict(formula="C9H10", density=1.03, fwhm_662=15.0, kind="scintillator"),
    "Si":        dict(formula="Si", density=2.33, fwhm_662=0.5, kind="semiconductor"),
    "HPGe":      dict(formula="Ge", density=5.32, fwhm_662=0.2, kind="semiconductor"),
    "CdTe":      dict(formula="CdTe", density=5.85, fwhm_662=2.0, kind="semiconductor"),
    "CZT":       dict(formula="Cd0.9Zn0.1Te", density=5.78, fwhm_662=2.0, kind="semiconductor"),
    "CsPbBr3":   dict(formula="CsPbBr3", density=4.55, fwhm_662=1.4, kind="semiconductor"),
    "FAPbBr3":   dict(formula="CH5N2PbBr3", density=3.8, fwhm_662=1.7, kind="semiconductor"),
}


def material(name, xcom, **overrides):
    """material('CZT', xcom, fwhm_662=1.5) – a crystal from the catalogue, with optional changes"""
    p = dict(CATALOGUE[name]); p.update(overrides)
    return Material(name=name, xcom=xcom, **p)


# --------------------------------------------------------------------------- sources
# gamma lines (keV) and emission probability per decay
SOURCES = {
    "Am-241": [(59.54, 0.359)],
    "Co-57":  [(122.06, 0.856), (136.47, 0.107), (14.41, 0.0916)],
    "Cs-137": [(661.66, 0.851)],
    "Co-60":  [(1173.2, 0.9985), (1332.5, 0.9998)],
    "Na-22":  [(511.0, 1.80), (1274.5, 0.999)],
    "K-40":   [(1460.8, 0.107)],
    "Tl-208": [(2614.5, 0.99), (583.2, 0.85)],
    "Ba-133": [(356.01, 0.62), (81.0, 0.33), (302.85, 0.18), (383.85, 0.089), (276.4, 0.072)],
}


# --------------------------------------------------------------------------- geometry
def _exit_dist(R, Lc, p, d):
    t = math.inf
    if d[0] > 1e-12: t = (Lc - p[0]) / d[0]
    elif d[0] < -1e-12: t = -p[0] / d[0]
    a = d[1] * d[1] + d[2] * d[2]
    if a > 1e-12:
        b = p[1] * d[1] + p[2] * d[2]; c = p[1] * p[1] + p[2] * p[2] - R * R
        tc = (-b + math.sqrt(max(0.0, b * b - a * c))) / a
        t = min(t, tc)
    return max(0.0, t)


def _entry_dist(R, Lc, p, d):
    c_ = []
    for xf in (0.0, Lc):
        if abs(d[0]) < 1e-12: continue
        t = (xf - p[0]) / d[0]
        if t <= 0: continue
        y, z = p[1] + d[1] * t, p[2] + d[2] * t
        if y * y + z * z <= R * R: c_.append(t)
    a = d[1] * d[1] + d[2] * d[2]
    if a > 1e-12:
        b = p[1] * d[1] + p[2] * d[2]; c = p[1] * p[1] + p[2] * p[2] - R * R
        disc = b * b - a * c
        if disc >= 0:
            for s in (-1, 1):
                t = (-b + s * math.sqrt(disc)) / a
                x = p[0] + d[0] * t
                if t > 0 and 0 <= x <= Lc: c_.append(t)
    return min(c_) if c_ else None


def _isotropic(rnd):
    c = 2 * rnd() - 1; s = math.sqrt(1 - c * c); p = 2 * math.pi * rnd()
    return (s * math.cos(p), s * math.sin(p), c)


def _rotate(d, cos, phi):
    sin = math.sqrt(max(0.0, 1 - cos * cos)); u, v, w = d
    if abs(w) > 0.99999:
        return (sin * math.cos(phi), sin * math.sin(phi), math.copysign(1, w) * cos)
    t = math.sqrt(1 - w * w)
    return (u * cos + sin * (u * w * math.cos(phi) - v * math.sin(phi)) / t,
            v * cos + sin * (v * w * math.cos(phi) + u * math.sin(phi)) / t,
            w * cos - sin * math.cos(phi) * t)


def _sample_compton(E, rnd):
    """Kahn's method: returns (E', cos θ)"""
    k = E / ME_C2
    while True:
        r1, r2, r3 = rnd(), rnd(), rnd()
        if r1 <= (1 + 2 * k) / (9 + 2 * k):
            x = 1 + 2 * k * r2
            if r3 <= 4 * (1 / x - 1 / (x * x)):
                return E / x, 1 - (x - 1) / k
        else:
            x = (1 + 2 * k) / (1 + 2 * k * r2)
            c = 1 - (x - 1) / k
            if r3 <= 0.5 * (c * c + 1 / x):
                return E / x, c


# --------------------------------------------------------------------------- one photon
@dataclass
class Setup:
    """Everything you can change about the measurement."""
    radius_cm: float = 3.15          # crystal radius (2.5" × 2.5" NaI: 3.15 cm)
    length_cm: float = 6.3           # crystal length along the axis (cm)
    distance_cm: float = 2.2         # source to front face (cm)
    window_mm: float = 1.0           # aluminium entrance window / housing (mm)
    shield: bool = True              # lead shield around the detector
    p_backscatter: float = 0.07      # shield: probability that a photon missing the crystal comes back backscattered
    p_lead_xray: float = 0.025       # shield: ... comes back as a Pb K X-ray (75 / 85 keV)
    p_annihilation: float = 0.025    # shield: ... comes back as a 511 keV annihilation photon (E > 1022 keV only)
    sampling: str = "demo"           # 'demo' = as on the website (most photons aimed at the crystal, 20 % isotropic
                                     #          with shield); 'cone' = isotropic emission inside the cone that
                                     #          contains the crystal (gives absolute efficiencies)


def run_photon(mat: Material, st: Setup, E0: float, rnd) -> tuple:
    """Follow one photon of energy E0. Returns (deposited energy keV, outcome).
    outcome: 'photo' (full energy), 'escape' (only a K X-ray left the crystal), 'compton', 'pair',
             'env' (came back from the shield), 'none'"""
    R, Lc = st.radius_cm, st.length_cm
    src = (-st.distance_cm, 0.0, 0.0)
    # initial direction
    if st.sampling == "cone":
        cos_max = st.distance_cm / math.hypot(st.distance_cm, R)
        c = 1 - rnd() * (1 - cos_max); s = math.sqrt(1 - c * c); ph = 2 * math.pi * rnd()
        d0 = (c, s * math.cos(ph), s * math.sin(ph))
    elif st.shield and rnd() < 0.2:
        d0 = _isotropic(rnd)
    else:   # aim at a random point of the front half of the crystal (website demo)
        rr = R * 1.15 * math.sqrt(rnd()); ph = 2 * math.pi * rnd()
        tgt = (Lc * 0.5 * rnd(), rr * math.cos(ph), rr * math.sin(ph))
        dd = (tgt[0] - src[0], tgt[1] - src[1], tgt[2] - src[2]); n = math.sqrt(sum(x * x for x in dd))
        d0 = (dd[0] / n, dd[1] / n, dd[2] / n)

    deposit = 0.0; any_pair = env = False; kx_emitted = []
    stack = [[src, d0, E0, False, False]]   # position, direction, energy, inside?, bounced?
    guard = 0
    while stack and guard < 60:
        guard += 1
        p, d, E, inside, bounced = stack.pop()
        if not inside:
            t = _entry_dist(R, Lc, p, d)
            if t is None:   # missed the crystal – maybe the shield sends something back
                if st.shield and not bounced and E > 40:
                    far = (p[0] + d[0] * 12, p[1] + d[1] * 12, p[2] + d[2] * 12)
                    u = rnd()
                    tgt = (Lc * rnd(), (rnd() - 0.5) * R, (rnd() - 0.5) * R)
                    b = (tgt[0] - far[0], tgt[1] - far[1], tgt[2] - far[2]); n = math.sqrt(sum(x * x for x in b))
                    nd = (b[0] / n, b[1] / n, b[2] / n)
                    if u < st.p_backscatter:
                        stack.append([far, nd, E / (1 + (E / ME_C2) * 1.95), False, True]); env = True
                    elif u < st.p_backscatter + st.p_lead_xray and E > 88:
                        stack.append([far, nd, 74.97 if rnd() < 0.8 else 84.9, False, True]); env = True
                    elif u < st.p_backscatter + st.p_lead_xray + st.p_annihilation and E > 1022:
                        stack.append([far, nd, 511.0, False, True]); env = True
                continue
            p = (p[0] + d[0] * t, p[1] + d[1] * t, p[2] + d[2] * t)
            mu_al = 2.7 * (0.075 + 3.4 * (20 / max(E, 10)) ** 2.9)
            if rnd() > math.exp(-mu_al * st.window_mm / 10):
                continue    # stopped in the housing
        for _ in range(40):  # inside the crystal
            mp, mc, mpp = mat.mu(E)
            mt = mp + mc + mpp
            s = -math.log(1 - rnd()) / mt
            te = _exit_dist(R, Lc, p, d)
            if s >= te:
                break        # escapes
            p = (p[0] + d[0] * s, p[1] + d[1] * s, p[2] + d[2] * s)
            u = rnd() * mt
            if u < mp:       # photo-absorption, maybe followed by a K X-ray
                kx = mat.fluorescence(E, rnd)
                if kx:
                    deposit += E - kx; kx_emitted.append(kx)
                    stack.append([p, _isotropic(rnd), kx, True, False])
                else:
                    deposit += E
                E = 0; break
            elif u < mp + mc:  # Compton
                E2, c = _sample_compton(E, rnd)
                deposit += E - E2; E = E2
                d = _rotate(d, c, 2 * math.pi * rnd())
                if E < 1:
                    deposit += E; E = 0; break
            else:            # pair production
                any_pair = True
                deposit += E - 1022
                a = _isotropic(rnd)
                stack.append([p, a, 511.0, True, False]); stack.append([p, (-a[0], -a[1], -a[2]), 511.0, True, False])
                E = 0; break

    if deposit <= 0:
        return 0.0, "none"
    full = abs(deposit - E0) < 0.5
    if env and not any_pair and deposit < E0 - 0.5: return deposit, "env"
    if full: return deposit, "photo"
    if any_pair: return deposit, "pair"
    if any(abs(deposit - (E0 - k)) < 0.5 for k in kx_emitted): return deposit, "escape"
    return deposit, "compton"


# --------------------------------------------------------------------------- a whole measurement
OUTCOMES = ("photo", "escape", "compton", "pair", "env")


def simulate(mat: Material, source, st: Setup | None = None, n_photons=10_000, emax=None, bins=400, seed=42):
    """Fire n_photons from `source` (a name in SOURCES, a list of (E, intensity) or one energy in keV).
    Returns a dict with the stacked histogram (counts per bin for each outcome) and counters."""
    st = st or Setup()
    lines = SOURCES[source] if isinstance(source, str) else ([(float(source), 1.0)] if isinstance(source, (int, float)) else source)
    tot_i = sum(i for _, i in lines)
    emax = emax or max(e for e, _ in lines) * 1.2 + 20
    rng = random.Random(seed); rnd = rng.random
    hist = {o: [0] * bins for o in OUTCOMES}
    interacted = peak = 0
    for _ in range(n_photons):
        u = rnd() * tot_i; acc = 0.0
        for E0, inten in lines:
            acc += inten
            if u <= acc: break
        dep, out = run_photon(mat, st, E0, rnd)
        if out == "none":
            continue
        interacted += 1
        if out == "photo": peak += 1
        e = dep + mat.fwhm(dep) / 2.355 * rng.gauss(0, 1)
        b = int(e / emax * bins)
        if 0 <= b < bins:
            hist[out][b] += 1
    return dict(material=mat.name, source=source if isinstance(source, str) else str(source), n=n_photons,
                interacted=interacted, peak=peak, emax=emax, bins=bins, hist=hist, setup=st.__dict__)
