import numpy as np
import math as m


def entrada_mp(
    t: float,
    RR: float,
    IEratio: float,
    Pmus_min: float,
    tau: float,
    Pao: float,
    Pvent: float = None,
) -> np.array:

    T = 60/RR
    Te = (60/RR)/(1 + IEratio)
    Ti = T - Te

    tciclo = t % T
    if tciclo <= Ti:  # inspiracao
        Pmus = (-Pmus_min / (Ti * Te)) * (tciclo**2) + ((Pmus_min * T) / (Ti * Te)) * tciclo
        dPmus = 2 * (-Pmus_min / (Ti * Te)) * tciclo + (Pmus_min * T / (Ti * Te))
    else:             # expiracao
        Pmus = (Pmus_min / (1 - m.exp(-Te / tau))) * (
            m.exp(-(tciclo - Ti) / tau) - m.exp(-Te / tau)
        )
        dPmus = (Pmus_min / (1 - m.exp(-Te / tau))) * ((-m.exp((Ti - tciclo) / tau)) / tau)

    u1 = dPmus
    u2 = Pvent if Pvent else Pao

    return np.array([u1, u2]), Pmus
