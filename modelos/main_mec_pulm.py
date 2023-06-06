from typing import List

import pandas as pd
import numpy as np
import math as m

from tqdm import tqdm

from parametros import cts_mp
from funcoes.derivada_mec_pulm import derivada_mp
from funcoes.entrada_mec_pulm import entrada_mp
from funcoes.saida_mec_pulm import saida_mp
from funcoes.plot_mec_pulm import plot_mp

from decorators.timefunc import timefunc

Pl = np.zeros(cts_mp["N"], dtype=int)
Ptr = np.zeros(cts_mp["N"], dtype=int)
Pb = np.zeros(cts_mp["N"], dtype=int)
PA = np.zeros(cts_mp["N"], dtype=int)
Ppl = np.zeros(cts_mp["N"], dtype=int)

Pao_Pvent_zero = np.zeros(cts_mp["N"], dtype=int)
dPmus_zero = np.zeros(cts_mp["N"], dtype=int)

vetor_zero = np.zeros(cts_mp["N"], dtype=int)

Ti = cts_mp["Ti"]
Te = cts_mp["Te"]
T = cts_mp["T"]
# Te = cts_mp["IEratio"]*Ti
tau = cts_mp["Te"]/5
RR = cts_mp["RR"]
IEratio = cts_mp["IEratio"]

Pmus_min = cts_mp["Pmus_min"]
Pao = cts_mp["Pao"]
# Pvent = cts_mp["Pvent"]
Pvent = None

Cl = cts_mp["Cl"]
Ctr = cts_mp["Ctr"]
Ccw = cts_mp["Ccw"]
Cb = cts_mp["Cb"]
CA = cts_mp["CA"]

Rtb = cts_mp["Rtb"]
Rlt = cts_mp["Rlt"]
Rml = cts_mp["Rml"]
RbA = cts_mp["RbA"]

Vul = cts_mp["Vul"]
Vut = cts_mp["Vut"]
Vub = cts_mp["Vub"]
VuA = cts_mp["VuA"]

f = RR/60
dt = cts_mp["dt"]


class MecanicaPulmonar:
    def __init__(self):
        self.t = np.arange(0, cts_mp["N"]*cts_mp["dt"], cts_mp["dt"])

        self.x = pd.DataFrame(
            {
                'Pl': Pl, 'Ptr': Ptr, 'Pb': Pb, 'PA': PA, 'Ppl': Ppl
            }
        )
        self.y = pd.DataFrame(
            {
                'V_dot': vetor_zero, 'VA_dot': vetor_zero, 'Vl': vetor_zero,
                'Vtr': vetor_zero, 'Vb': vetor_zero, 'VA': vetor_zero,
                'VD': vetor_zero, 'V': vetor_zero,
            }
        )
        self.u = pd.DataFrame(
            {
                'dPmus': vetor_zero, 'Pao_Pvent': vetor_zero, 'Pmus': vetor_zero
            }
        )

    @timefunc
    def run_mecanica_pulmonar(self):
        """
            The initial conditions for the five different
            pressure nodes in Fig. 6 (state variables) have been assigned, assuming that at time t=0,
            corresponding to the end-exhalation time, all of the pressures in the lungs equilibrate to Patm,
            whereas the Ppl has a subatmospheric value of -5 cmH2O (72)
        """
        self.x.iloc[0, 0] = Pao
        self.x.iloc[0, 1] = Pao
        self.x.iloc[0, 2] = Pao
        self.x.iloc[0, 3] = Pao
        self.x.iloc[0, 4] = -5
        self.rungekutta4()

    def plot_mecanica_pulmonar(self):
        plot_mp(self.t, self.x, self.y, self.u)

    def rungekutta4(self):
        for i in tqdm(range(cts_mp["N"]-1)):
            # u em t
            t = self.t[i]
            # T, Ti, Te, RR = controle_mp(t, RR, PO2, PCO2...)

            phi = 2 * m.pi * f * t
            u_t_array, Pmus = entrada_mp(t, RR, IEratio, Pmus_min, tau, Pao, Pvent)
            self.u.loc[i, 'dPmus'] = u_t_array[0]
            self.u.loc[i, 'Pao_Pvent'] = u_t_array[1]
            self.u.loc[i, 'Pmus'] = Pmus
            u_t = np.matrix(u_t_array).transpose()

            # u em t+dt/2
            t_meio = (t+(t+dt))/2
            phi_meio = 2 * m.pi * f * t_meio
            u_meio_array, Pmus_meio = entrada_mp(t_meio, RR, IEratio, Pmus_min, tau, Pao, Pvent)
            u_meio = np.matrix(u_meio_array).transpose()

            # u em t+dt*2/3
            t_23 = (t+dt*(2/3))
            phi_23 = 2 * m.pi * f * t_23
            u_23_array, Pmus_23 = entrada_mp(t_23, RR, IEratio, Pmus_min, tau, Pao, Pvent)
            u_23 = np.matrix(u_23_array).transpose()

            # u em t+dt
            t_dt = t + dt
            phi_dt = 2 * m.pi * f * t_dt
            u_dt_array, Pmus_dt = entrada_mp(t_dt, RR, IEratio, Pmus_min, tau, Pao, Pvent)
            u_dt = np.matrix(u_dt_array).transpose()

            # x em t
            x_array = self.x.iloc[i, 0:].to_numpy()
            x = np.matrix(x_array).transpose()

            # constantes para RK
            k1 = derivada_mp(x, u_t,  Cl, Ctr, Ccw, Cb, CA, Rtb, Rlt, Rml, RbA)
            k2 = derivada_mp(x + k1*dt/2, u_meio, Cl, Ctr, Ccw, Cb, CA, Rtb, Rlt, Rml, RbA)
            k3 = derivada_mp(x + k2*dt/2, u_meio, Cl, Ctr, Ccw, Cb, CA, Rtb, Rlt, Rml, RbA)
            k4 = derivada_mp(x + k3*dt, u_dt, Cl, Ctr, Ccw, Cb, CA, Rtb, Rlt, Rml, RbA)

            # constantes para RK2
            k1_rk2 = derivada_mp(x, u_t, Cl, Ctr, Ccw, Cb, CA, Rtb, Rlt, Rml, RbA)
            k2_rk2 = derivada_mp(x + k1_rk2*dt, u_dt, Cl, Ctr, Ccw, Cb, CA, Rtb, Rlt, Rml, RbA)

            # atribuindo os valores em x
            # jogar em variaveis intermediarias antes de transpor
            # x = x + dt*(k1+2*k2+2*k3+k4)/6
            # x_rk2 = x + dt*(k1_rk2*0.25+k2_rk2*0.75)
            x_rk2 = x + dt*(k1_rk2+k2_rk2)/2

            Pl = x[0, 0]
            Ptr = x[1, 0]
            Pb = x[2, 0]
            PA = x[3, 0]
            Ppl = x[4, 0]
            # self.x.iloc[i+1, :] = pd.DataFrame(x.transpose())
            self.x.iloc[i+1, :] = pd.DataFrame(x_rk2.transpose())
            # self.x.iloc[i+1, 0] = x[0, 0] + dt*(k1[0, 0]+2*k2[0, 0]+2*k3[0, 0]+k4[0, 0])/6
            # self.x.iloc[i+1, 1] = x[1, 0] + dt*(k1[1, 0]+2*k2[1, 0]+2*k3[1, 0]+k4[1, 0])/6
            # self.x.iloc[i+1, 2] = x[2, 0] + dt*(k1[2, 0]+2*k2[2, 0]+2*k3[2, 0]+k4[2, 0])/6
            # self.x.iloc[i+1, 3] = x[3, 0] + dt*(k1[3, 0]+2*k2[3, 0]+2*k3[3, 0]+k4[3, 0])/6
            # self.x.iloc[i+1, 4] = x[4, 0] + dt*(k1[4, 0]+2*k2[4, 0]+2*k3[4, 0]+k4[4, 0])/6

            # atribuindo os valores em y
            self.y.iloc[i+1, :] = saida_mp(x, u_t, Cl, Ctr, Cb, CA, Rtb, Rlt, Rml, RbA, Vul, Vut, Vub, VuA)

            phi = phi_dt
