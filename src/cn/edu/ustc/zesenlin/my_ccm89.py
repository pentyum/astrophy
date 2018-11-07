import pylab as pl

def get_Fa_Fb(x):
    if (x >= 5.9) and (x <= 8.):
        Fa = -0.04473 * (x - 5.9)**2 - 0.009779 * (x - 5.9)**3
        Fb = 0.2130 * (x - 5.9)**2 + 0.1207 * (x - 5.9)**3
    else:
        Fa = 0.
        Fb = 0.
    return Fa, Fb


def get_a_b(x):
    if (x >= 0.3) and (x <= 1.1):
        a = 0.574 * x**1.61
        b = -0.527 * x**1.61
    elif (x > 1.1) and (x <= 3.3):
        y = x - 1.82
        a = 1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 0.72085 * y**4 +\
            0.01979 * y**5 - 0.77530 * y**6 + 0.32999 * y**7
        b = 1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 5.38434 * y**4 -\
            0.62251 * y**5 + 5.30260 * y**6 - 2.09002 * y**7
    elif (x > 3.3) and (x <= 8.):
        Fa, Fb = get_Fa_Fb(x)
        a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67)**2 + 0.341) + Fa
        b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62)**2 + 0.263) + Fb
    elif (x > 8.) and (x < 10.):
        a = -1.073 - 0.628 * (x - 8.) + 0.137 * (x - 8.)**2 - 0.070 * (x - 8.)**3
        b = 13.670 + 4.257 * (x - 8.) - 0.420 * (x - 8.)**2 + 0.374 * (x - 8.)**3
    else:
        a = pl.nan
        b = pl.nan
    return a, b


def ccm89(wav, a_v, r_v=3.1):
    '''Get dust extincetion using CCM89 law.

    Author: ZSLIN (zesenlin@mail.ustc.edu.cn)'''
    x = 1. / wav * 1.e4  # convert Angstrom to 1./micron
    nwav = len(wav)
    A_lam = pl.zeros_like(wav)
    for i in range(nwav):
        a, b = get_a_b(x[i])
        A_lam[i] = (a + b / r_v) * a_v
    return A_lam