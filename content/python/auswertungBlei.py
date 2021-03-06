import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.constants as const
from uncertainties import ufloat

cw = 4.18

roh, M, alpha, kappa = np.genfromtxt('../values/material_val.txt', unpack=True)

roh=roh[2]
M=M[2]
alpha=alpha[2]
kappa=kappa[2]

roh *= 1e6
alpha *= 1e-6
kappa *= 1e6

aml, amg, amwl, aTk, aTw, aTg = np.genfromtxt('../values/standard_val.txt', unpack=True)
amw = amwl - aml
amk = amg - amwl

x_, x_, mal, x_, x_, nmal = np.genfromtxt('../values/mass_val.txt', unpack=True)
mass = mal - nmal

Tw, Ts, Te, mw = np.genfromtxt('../values/blei_val.txt', unpack=True)
mTw = np.mean(Tw)
sTw = np.std(Tw)
fTw = ufloat(mTw, sTw)
mTs = np.mean(Ts)
sTs = np.std(Ts)
fTs = ufloat(mTs, sTs)
mTe = np.mean(Te)
sTe = np.std(Te)
fTe = ufloat(mTe, sTe)
mmw = np.mean(mw)
smw = np.std(mw)
fmw = ufloat(mmw, smw)

aTk = const.convert_temperature(aTk,'c','K')
aTw = const.convert_temperature(aTw,'c','K')
aTg = const.convert_temperature(aTg,'c','K')
fTw = const.convert_temperature(fTw,'c','K')
fTs = const.convert_temperature(fTs,'c','K')
fTe = const.convert_temperature(fTe,'c','K')

def cgmg(Tx, Ty, Tm, cw, mx, my):
    return (cw*my*(Ty-Tm)-cw*mx*(Tm-Tx))/(Tm-Tx)

def ck(cgmg, cw, mw, mk, Tm, Tw, Tk):
    return ((cw*mw+cgmg)*(Tm-Tw))/(mk*(Tk-Tm))

def atomwaerme(alpha, kappa, M, roh, Tm, ck):
    return (ck * M) - (9 * (alpha**2) * kappa * (M / roh) * Tm)

ck = ck(cgmg(aTk, aTw, aTg, cw, amk, amw), cw, fmw, mass, fTe, fTw, fTs)
CV = atomwaerme(alpha, kappa, M, roh, fTe, ck)

print('Blei: ck =',ck,'Joule/(g*K)     CV = ',CV,'Joule/(mol*K)       Abweichung: ',(3*8.314 - CV) / (3*8.314) * 100,'%')
