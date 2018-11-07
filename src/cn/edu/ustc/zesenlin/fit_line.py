import pylab as pl
from astropy.io import ascii
from cn.edu.ustc.zesenlin.my_ccm89 import ccm89
from com.google.code.astrolibpy.mpfit.mpfit3 import mpfit


def deredden(wav, flux, err=None, ebv=0.):
    R_V = 3.1
    A_V = ebv * R_V
    A_lam = ccm89(wav, A_V, R_V)
    flux_deredden = flux * 10**(-0.4 * A_lam)
    if err is None:
        return flux_deredden
    else:
        err_deredden = err * 10**(-0.4 * A_lam)
        return flux_deredden, err_deredden


def gaussian(x, mu, sigma, A):
    '''Gaussian function.'''
    y = A * pl.exp(-(x - mu)**2 / (2. * sigma**2)) / pl.sqrt(2. * pl.pi * sigma**2)
    return y


def single_gaussian(x, p):
    '''Function form of single gaussian profile.'''
    y = p[0] + gaussian(x, *p[1:])
    return y


def double_gaussian(x, p):
    '''Function form of double gaussian profiles.'''
    y = p[0] + gaussian(x, *p[1:4]) + gaussian(x, *p[4:7])
    return y


def tri_gaussian(x, p):
    '''Function form of three gaussian profiles.'''
    y = p[0] + gaussian(x, *p[1:4]) + gaussian(x, *p[4:7]) + gaussian(x, *p[7:10])
    return y


def myfunct_gaussian(p, fjac=None, x=None, y=None, err=None, ngaussian=1):
    # Parameter values are passed in "p"
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.
    if ngaussian == 1:
        model = single_gaussian(x, p)
    elif ngaussian == 2:
        model = double_gaussian(x, p)
    elif ngaussian == 3:
        model = tri_gaussian(x, p)
    else:
        print('nguassian must be smaller than 3.')
        return [-1, pl.NaN]
    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    status = 0
    return [status, (y - model) / err]


def generate_parinfo(p0, plimits=None, pname=None, ltype='emi'):
    '''Generate the parameter information in mpfit format for each line.'''
    parinfo = [{'value': 0., 'fixed': 0, 'limited': [1, 1], 'limits': [0., 0.], 'tied': '', 'parname': 'name'} for i in range(3)]
    if plimits is None:
        plimits = [5., 1., 0.]
    if pname is None:
        pname = ['wav_center', 'width', 'area']
    for i in range(3):
        parinfo[i]['value'] = p0[i]
        parinfo[i]['limits'] = [p0[i] - plimits[i], p0[i] + plimits[i]]
        parinfo[i]['parname'] = pname[i]
    if ltype == 'emi':
        parinfo[2]['limits'] = [0., pl.Inf]
    elif ltype == 'abs':
        parinfo[2]['limits'] = [-pl.Inf, 0.]
    return parinfo


class LineProfile:
    def __init__(self, name, wavelength, flux, flux_err, sigma, sigma_err, base):
        self.name = name
        self.wavelength = wavelength
        self.flux = flux
        self.flux_err = flux_err
        self.sigma = sigma
        self.sigma_err = sigma_err
        self.snr = None
        self.ew = None
        self.base = base
        self.spec = None

    def profile(self):
        if self.spec is None:
            print('Warning: original spectrum is not saved!')
        else:
            pl.plot(self.spec[0], self.spec[1], 'k-')
            pl.plot(self.spec[0], single_gaussian(self.spec[0], [self.base, self.wavelength, self.sigma, self.flux]), 'r-')
            pl.show()


def emptyline(lname, wav_center):
    try:
        nline = len(wav_center)
    except TypeError:
        empty = LineProfile(lname, wav_center, 0., 0., 0., 0., 0.)
    else:
        empty = [LineProfile(lname[i], wav_center[i], 0., 0., 0., 0., 0.) for i in range(nline)]
    return empty


def fit_single_gaussian(wav, flux, flux_err, name_line, wav_line, isplot=False, silent=True):
    p0 = [pl.median(flux), wav_line, 2., 10.]
    parinfo = [{'value': 0., 'fixed': 0, 'limited': [0, 0], 'limits': [0., 0.], 'parname': 'name'}]
    parinfo[0]['value'] = p0[0]
    parinfo[0]['parname'] = 'constant'
    parinfo.extend(generate_parinfo(p0=p0[1:]))
    fdata = {'x': wav, 'y': flux, 'err': flux_err, 'ngaussian': 1}
    res = mpfit(myfunct_gaussian, p0, parinfo=parinfo, functkw=fdata, quiet=silent)
    if (res.status < 0) or (res.perror is None):
        print('error message = ', res.errmsg)
        return emptyline(name_line, wav_line)
    line = LineProfile(name_line, res.params[1], res.params[3], res.perror[3], res.params[2], res.perror[2], res.params[0])
    if isplot:
        line.spec = [wav, flux]
    return line


def fit_double_gaussian(wav, flux, flux_err, name_line, wav_line, rline=None, isplot=False, silent=True):
    p0 = [pl.median(flux), wav_line[0], 2., 10., wav_line[1], 2., 10.]
    parinfo = [{'value': 0., 'fixed': 0, 'limited': [0, 0], 'limits': [0., 0.], 'tied': '', 'parname': 'name'}]
    parinfo[0]['value'] = p0[0]
    parinfo[0]['parname'] = 'constant'
    parinfo.extend(generate_parinfo(p0=p0[1:4], pname=['wav_center1', 'width1', 'area1']))
    parinfo.extend(generate_parinfo(p0=p0[4:7], pname=['wav_center2', 'width2', 'area2']))
    if rline is not None:
        parinfo[3]['tied'] = str(rline) + ' * p[6]'
    fdata = {'x': wav, 'y': flux, 'err': flux_err, 'ngaussian': 2}
    res = mpfit(myfunct_gaussian, p0, parinfo=parinfo, functkw=fdata, quiet=silent)
    if (res.status < 0) or (res.perror is None):
        print('error message = ', res.errmsg)
        return emptyline(name_line, wav_line)
    line1 = LineProfile(name_line[0], res.params[1], res.params[3], res.perror[3], res.params[2], res.perror[2], res.params[0])
    line2 = LineProfile(name_line[1], res.params[4], res.params[6], res.perror[6], res.params[5], res.perror[5], res.params[0])
    if isplot:
        line1.spec = [wav, flux]
        line2.spec = [wav, flux]
    return [line1, line2]


def fit_tri_gaussian(wav, flux, flux_err, name_line, wav_line, rline=None, isplot=False, silent=True):
    p0 = [pl.median(flux), wav_line[0], 2., 10., wav_line[1], 2., 10., wav_line[2], 2., 10.]
    parinfo = [{'value': 0., 'fixed': 0, 'limited': [0, 0], 'limits': [0., 0.], 'tied': '', 'parname': 'name'}]
    parinfo[0]['value'] = p0[0]
    parinfo[0]['parname'] = 'constant'
    parinfo.extend(generate_parinfo(p0=p0[1:4], pname=['wav_center1', 'width1', 'area1']))
    parinfo.extend(generate_parinfo(p0=p0[4:7], pname=['wav_center2', 'width2', 'area2']))
    parinfo.extend(generate_parinfo(p0=p0[7:10], pname=['wav_center3', 'width3', 'area3']))
    if rline is not None:
        if rline[0] is not None:
            parinfo[3]['tied'] = str(rline[0]) + ' * p[9]'
        if rline[1] is not None:
            parinfo[6]['tied'] = str(rline[1]) + ' * p[9]'
    fdata = {'x': wav, 'y': flux, 'err': flux_err, 'ngaussian': 3}
    res = mpfit(myfunct_gaussian, p0, parinfo=parinfo, functkw=fdata, quiet=silent)
    if (res.status < 0) or (res.perror is None):
        print('error message = ', res.errmsg)
        return emptyline(name_line, wav_line)
    line1 = LineProfile(name_line[0], res.params[1], res.params[3], res.perror[3], res.params[2], res.perror[2], res.params[0])
    line2 = LineProfile(name_line[1], res.params[4], res.params[6], res.perror[6], res.params[5], res.perror[5], res.params[0])
    line3 = LineProfile(name_line[2], res.params[7], res.params[9], res.perror[9], res.params[8], res.perror[8], res.params[0])
    if isplot:
        line1.spec = [wav, flux]
        line2.spec = [wav, flux]
        line3.spec = [wav, flux]
    return [line1, line2, line3]


def multi_line_info(line_id):
    line_info = {
        17: [[17, 18], None],
        18: [[17, 18], None],
        29: [[29, 30], 1. / 3.],
        30: [[29, 30], 1. / 3.],
        45: [[45, 46, 47], [0.348116, None]],
        46: [[45, 46, 47], [0.348116, None]],
        47: [[45, 46, 47], [0.348116, None]],
        49: [[49, 50], None],
        50: [[49, 50], None]
    }
    return line_info.get(line_id, None)


def fit_line(wav, flux, flux_err=None, lines=[0, 13, 17, 25, 29, 45, 49], path_line='./', isplot=False, silent=True):
    '''Fit a given emission/absorption line set for a rest-frame spectrum.
    Input
    --------
    wav: array like, wavelength, in units of Angstrom, in rest-frame.
    flux: flux array, in rest-frame.
    flux_err: error of flux, in rest-frame. If flux_err is not given, error will be set to 0.1*flux.
    lines: The index set (in lines.dat file) that which lines to be fitted, array like.
    path_line: str, path of lines.dat file.
    isplot: whethere or not to plot the fitting result. If True, a method named profile will be saved in the return line class.

    Return
    --------
    res: list of LineProfile class (see fit_line.py).
         Each element gives properties (i.e., name, wavelength, flux, flux error, sigma, sigma error,
         and base in fitting) of one line. The value of these properties can be extracted like:
         res[0].name, res[0].wavelength, res[0].flux, res[0].flux_err, res[0].sigma, res[0].sigma_err,
         res[0].base.

    Author: ZSLIN (zesenlin@mail.ustc.edu.cn)'''

    # if flux error is not given set to 10% of flux
    if flux_err is None:
        flux_err = flux * 0.1
    # read the line table
    fline = ascii.read(path_line + 'lines.dat')
    line_wav = fline['Wavelength']
    line_name = fline['Line']
    # only use spectrum piece with width of sp_bin*2 to fit
    sp_bin = 25.  # in AA
    nline = len(lines)
    res_lines = []
    for i in range(nline):
        line_info = multi_line_info(lines[i])
        if line_info is None:
            wav_center = line_wav[lines[i]]
            lname = line_name[lines[i]]
            isin = abs(wav - wav_center) <= sp_bin
            wav_use = wav[isin]
            if (len(wav_use) > 0) and (wav_center > min(wav_use)) and (wav_center < max(wav_use)):
                res_fit = fit_single_gaussian(wav_use, flux[isin], flux_err[isin], name_line=lname, wav_line=wav_center, isplot=isplot)
            else:
                res_fit = emptyline(lname, wav_center)
            res_lines.append(res_fit)
        elif len(line_info[0]) == 2:
            wav_center = [line_wav[ind] for ind in line_info[0]]
            line_ratio = line_info[1]
            isin = (wav >= wav_center[0] - sp_bin) * (wav <= wav_center[-1] + sp_bin)
            lname = [line_name[ind] for ind in line_info[0]]
            if (len(wav[isin]) > 0):
                res_fit = fit_double_gaussian(wav[isin], flux[isin], flux_err[isin], name_line=lname, wav_line=wav_center,
                                              rline=line_ratio, isplot=isplot)
            else:
                res_fit = emptyline(lname, wav_center)
            res_lines.extend(res_fit)
        elif len(line_info[0]) == 3:
            wav_center = [line_wav[ind] for ind in line_info[0]]
            line_ratio = line_info[1]
            isin = (wav >= wav_center[0] - sp_bin) * (wav <= wav_center[-1] + sp_bin)
            lname = [line_name[ind] for ind in line_info[0]]
            if (len(wav[isin]) > 0):
                res_fit = fit_tri_gaussian(wav[isin], flux[isin], flux_err[isin], name_line=lname, wav_line=wav_center, rline=line_ratio,
                                           isplot=isplot)
            else:
                res_fit = emptyline(lname, wav_center)
            res_lines.extend(res_fit)
    return res_lines