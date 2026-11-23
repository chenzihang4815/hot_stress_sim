# -*- coding: utf-8 -*-
"""Temperature-dependent conductivity and Kirchhoff transform.

kappa(T) = kappa0 * (T_ref / T)^alpha

Silicon in the paper uses alpha ~ 1.5; copper uses alpha ~ 0.
The Kirchhoff variable is defined per the paper:

    psi = 1/kappa(T_g) * integral_{T_g}^{T} kappa(eta) d eta

so that kappa(T) grad T = kappa(T_g) grad psi inside one Picard step.
"""

from __future__ import print_function

import numpy as np


class Material(object):
    def __init__(self, name, kappa0, T_ref=300.0, alpha=0.0, rho=2330.0, cp=700.0,
                 E=None, nu=0.3, alpha_cte=0.0):
        self.name = name
        self.kappa0 = float(kappa0)
        self.T_ref = float(T_ref)
        self.alpha = float(alpha)
        self.rho = float(rho)
        self.cp = float(cp)
        self.E = float(E) if E is not None else None
        self.nu = float(nu)
        self.alpha_cte = float(alpha_cte)

    @property
    def volumetric_heat_capacity(self):
        return self.rho * self.cp

    def kappa(self, T):
        T = np.asarray(T, dtype=float)
        T_safe = np.maximum(T, 50.0)
        if self.alpha == 0.0:
            if T_safe.shape == ():
                return self.kappa0
            return np.full(T_safe.shape, self.kappa0, dtype=float)
        return self.kappa0 * (self.T_ref / T_safe) ** self.alpha

    def kirchhoff_forward(self, T, T_g):
        """Map physical temperature T to psi, referenced at T_g."""
        T = np.asarray(T, dtype=float)
        T_g = np.asarray(T_g, dtype=float)
        T = np.maximum(T, 50.0)
        T_g = np.maximum(T_g, 50.0)
        a = self.alpha
        if abs(a) < 1e-12:
            return T - T_g
        if abs(a - 1.0) < 1e-12:
            return T_g * np.log(T / T_g)
        return (T_g ** a) / (1.0 - a) * (T ** (1.0 - a) - T_g ** (1.0 - a))

    def kirchhoff_inverse(self, psi, T_g):
        """Recover T from psi. Vectorized."""
        psi = np.asarray(psi, dtype=float)
        T_g = np.asarray(T_g, dtype=float)
        T_g = np.maximum(T_g, 50.0)
        a = self.alpha
        if abs(a) < 1e-12:
            return T_g + psi
        if abs(a - 1.0) < 1e-12:
            return T_g * np.exp(psi / T_g)
        inner = T_g ** (1.0 - a) + psi * (1.0 - a) / (T_g ** a)
        inner = np.maximum(inner, 1e-12)
        T = inner ** (1.0 / (1.0 - a))
        return np.clip(T, 50.0, 2000.0)


def material_from_dict(name, d):
    return Material(
        name=name,
        kappa0=d["kappa0_W_mK"],
        T_ref=d.get("T_ref_K", 300.0),
        alpha=d.get("alpha", 0.0),
        rho=d.get("rho_kg_m3", 2330.0),
        cp=d.get("cp_J_kgK", 700.0),
        E=d.get("E_Pa"),
        nu=d.get("nu", 0.3),
        alpha_cte=d.get("alpha_cte_1_K", 0.0),
    )


def harmonic_mean(k1, k2):
    s = k1 + k2
    return np.where(s > 0.0, 2.0 * k1 * k2 / s, 0.0)
