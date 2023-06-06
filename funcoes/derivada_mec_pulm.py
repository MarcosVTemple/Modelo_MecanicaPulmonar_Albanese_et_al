import numpy as np


def derivada_mp(
    x: np.ndarray,
    u: np.ndarray,
    Cl: float,
    Ctr: float,
    Ccw: float,
    Cb: float,
    CA: float,
    Rtb: float,
    Rlt: float,
    Rml: float,
    RbA: float,
):
    Pl, Ptr, Pb, PA, Ppl = x[:, 0]
    dPmus, Pao = u[:, 0]
    dPl = (1 / (Cl * Rml))*Pao + (-(1 / (Cl * Rml)) - (1 / (Cl * Rlt)))*Pl + (1 / (Cl * Rlt))*Ptr
    dPtr = ((1 / (Ctr * Rlt)) + (1 / (Ccw * Rlt)))*Pl + (-(1 / (Ctr * Rlt)) - (1 / (Ctr * Rtb)) - (1 / (Ccw * Rlt)))*Ptr + (1 / (Ctr * Rtb))*Pb + dPmus
    dPb = (1 / (Ccw * Rlt))*Pl + ((1 / (Cb * Rtb)) - (1 / (Ccw * Rlt)))*Ptr + (-(1 / (Cb * Rtb)) - (1 / (Cb * RbA)))*Pb + (1 / (Cb * RbA))*PA + dPmus
    dPA = (1 / (Ccw * Rlt))*Pl + (-(1 / (Ccw * Rlt)))*Ptr + (1 / (CA * RbA))*Pb + (-(1 / (CA * RbA)))*PA + dPmus
    dPpl = (1 / (Ccw * Rlt))*Pl + (-(1 / (Ccw * Rlt)))*Ptr + dPmus

    a = np.matrix(
        [
            [
                -(1 / (Cl * Rml)) - (1 / (Cl * Rlt)), (1 / (Cl * Rlt)), 0, 0, 0
            ],
            [
                (1 / (Ctr * Rlt)) + (1 / (Ccw * Rlt)), -(1 / (Ctr * Rlt)) - (1 / (Ctr * Rtb)) - (1 / (Ccw * Rlt)), (1 / (Ctr * Rtb)), 0, 0,
            ],
            [
                (1 / (Ccw * Rlt)), (1 / (Cb * Rtb)) - (1 / (Ccw * Rlt)), -(1 / (Cb * Rtb)) - (1 / (Cb * RbA)), (1 / (Cb * RbA)), 0,
            ],
            [
                (1 / (Ccw * Rlt)), -(1 / (Ccw * Rlt)), (1 / (CA * RbA)), -(1 / (CA * RbA)), 0,
            ],
            [
                (1 / (Ccw * Rlt)), -(1 / (Ccw * Rlt)), 0, 0, 0
            ],
        ],
    )
    b = np.matrix(
        [
            [0, (1 / (Cl * Rml))],
            [1, 0],
            [1, 0],
            [1, 0],
            [1, 0]
        ]
    )

    # dx = np.dot(a, x) + np.dot(b, u) #8min19s
    dx = np.matrix(np.array([dPl, dPtr, dPb, dPA, dPpl])).transpose() #8min19s
    return dx
