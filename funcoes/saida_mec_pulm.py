import numpy as np
from numpy import ndarray


def saida_mp(
    x: np.ndarray,
    u: np.ndarray,
    Cl: float,
    Ctr: float,
    Cb: float,
    CA: float,
    Rtb: float,
    Rlt: float,
    Rml: float,
    RbA: float,
    Vul: float,
    Vut: float,
    Vub: float,
    VuA: float,
) -> ndarray:
    Pl, Ptr, Pb, PA, Ppl = x[:, 0]
    dPmus, Pao = u[:, 0]

    V_dot = -(1 / Rml) * Pl + (1 / Rml) * Pao
    VA_dot = (1 / RbA) * Pb + (-1 / RbA) * PA
    Vl = Cl * Pl + Vul
    Vtr = Ctr * Ptr - Ctr * Ppl + Vut
    Vb = Cb * Pb - Cb * Ppl + Vub
    VA = CA * PA - CA * Ppl + VuA
    VD = Vl + Vtr + Vb
    V = VD + VA

    y = np.array([V_dot, VA_dot, Vl, Vtr, Vb, VA, VD, V])
    # y = np.matrix(y_array).transpose()
    return y
