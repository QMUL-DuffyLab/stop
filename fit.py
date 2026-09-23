# -*- coding: utf-8 -*-
import os            # For file and directory handling
import numpy as np    # For numerical operations
import pandas as pd   # For data manipulation with DataFrames
import matplotlib.pyplot as plt  # For plotting
from scipy.optimize import nnls, differential_evolution

'''
fitting of TCSPC curves output from STOP using
NNLS and differential evolution. code mostly taken from
https://github.com/mkizilov/ReconFit
with unneeded stuff removed
'''

def gen_irf(filename, x):
    '''
    use the input file for the fortran to pull a FWHM for the pulse (IRF)
    and generate the corresponding curve. in future, add this to the
    fortran histogram!
    '''
    with open(filename, "r") as f:
        lines = [line.rstrip() for line in f]
    fwhm = float(lines[0])
    sigma = fwhm / 2.355
    mu = fwhm * 2.0 # fortran puts the peak there - will get shifted anyway
    irf = ((1 / (sigma * np.sqrt(2. * np.pi))) *
            np.exp(-(x - mu)**2 / (np.sqrt(2.) * sigma)**2))
    irf /= np.sum(irf)
    return irf

def get_histogram(filename, irf_file=None):
    '''
    take the output from the fortran and return the stuff
    we need to plot and fit it all. pandas/xarray do not like this
    because there are basically two header rows which we need for
    different things - the labels for each decay pathway and then
    the row telling us which are emissive. so just do it manually
    '''
    with open(filename, "r") as f:
        lines = [line.rstrip() for line in f]
    labels = lines[0].split(" ")
    labels.append("Emitted")
    emissive_str = lines[1].split(" ")
    emissive = [True if s == "T" else False for s in emissive_str]
    str_array = [line.split(" ") for line in lines[2:]]
    data = np.array([row for row in str_array]).astype(float)
    x = data[:, 0]
    counts = data[:, 1:]
    emissive_counts = counts[:, emissive[1:]]
    if emissive_counts.shape[1] > 1:
        sum_emissive = np.sum(emissive_counts, axis=0)
    else:
        sum_emissive = emissive_counts
    all_columns = np.hstack((data, sum_emissive))
    d = {l: c for l, c in zip(labels, all_columns.T)}
    if 'IRF' not in labels:
        sim_file = os.path.join(
                os.path.dirname(filename), "simulation_params")
        d['IRF'] = gen_irf(sim_file, d['Time(s)'])
    df = pd.DataFrame(d)
    return df

def get_si_exponent(x, latex = True):
    '''
    return the SI exponent of a number and short/long prefixes
    '''
    # "μ" doesn't work with latex
    if latex:
        mustr = r'$ \mu $'
    else:
        mustr = "μ"
    exponent = np.floor(np.floor(np.log10(float(x))) / 3) * 3
    pref = {-12: ["p", "pico"], -9: ["n", "nano"],
            -6: [mustr, "micro"], -3: ["m", "milli"], 0: ["", ""],
            3: ["k", "kilo"], 6: ["M", "mega"], 9: ["G", "giga"],
            12: ["T", "tera"]}
    if exponent in pref:
        return (exponent, pref[exponent])
    else:
        return (exponent, [None, None])

def sci_format(x):
    exponent, prefs = get_si_exponent(x)
    if exponent == 0:
        xp = f"{x / 10**exponent:.4f} " + "s"
    else:
        xp = f"{x / 10**exponent:.4f} " + prefs[0] + "s"
    return xp

def fix_x(x, nticks, ax):
    '''
    sort out the x axis of a given matplotlib.axis object
    '''
    xmax = np.max(x)
    xticks = [i * (xmax) / (nticks - 1) for i in range(nticks)]
    exponent, prefs = get_si_exponent(xmax)
    xlabel = "Time (" + prefs[0] + "s)"
    xtickround = np.around([i * (xmax / 10**exponent) / (nticks - 1)
        for i in range(nticks)]).astype(int)
    xticklabels = [str(i) for i in xtickround]
    ax.set_xlabel(xlabel)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    return ax

def fix_y(ax, ymax):
    '''
    sort out the y axis of a given matplotlib.axis object
    '''
    ax.set_yscale('log')
    ax.set_ylim([1e-5 * ymax, 1.1 * ymax])
    if ymax == 1.:
        ax.set_ylabel("Counts (normalised)")
    else:
        ax.set_ylabel("Counts")
    return ax

def plot_setup(x, nticks, ymax, cfigax=None):
    '''
    set up the figure and axes since they're gonna be the
    same for all the plots we do here
    '''
    if cfigax is None:
        fig, ax = plt.subplots(figsize=(12,8))
        plt.grid(visible=True)
    else:
        fig, ax = cfigax
    ax = fix_x(x, nticks, ax)
    ax = fix_y(ax, ymax)
    return fig, ax
    
def plot_all(df, outfile):
    '''
    big plot with all binned decay pathways plotted, including
    ones which would be invisible to a real detector.
    '''
    fig, ax = plot_setup(df['Time(s)'], 4, np.max(df['Emitted']))
    for col in df.columns:
        if np.sum(df[col]) > 0:
            plt.plot(df['Time(s)'], df[col], label=col)
    ax.legend()
    ax.set_ylim([0.01, 1.1 * np.max(df.values)])
    fig.tight_layout()
    plt.savefig(outfile)
    plt.close()

def plot_all_from_file(histfile):
    '''
    wrap get_histogram and plot_all to make it easier
    if plotting bits in a terminal or whatever
    '''
    df = get_histogram(histfile)
    outfile = os.path.splitext(histfile)[0] + ".pdf"
    plot_all(df, outfile)

def convolve(x, h):
    X = np.fft.fft(x, n=len(x) + len(h) - 1)
    H = np.fft.fft(h, n=len(x) + len(h) - 1)
    xch = np.real(np.fft.ifft(X * H))
    return xch[:len(x)]

def shift_irf(irf, irf_shift):
    n = len(irf)
    channel = np.arange(n)
    frac = irf_shift - np.floor(irf_shift)
    back = np.fmod(channel - np.floor(irf_shift), n)
    fwd  = np.fmod(channel - np.ceil(irf_shift), n)
    # contribution from lower bin
    broll = (1 - frac) * irf[np.fmod(back + n, n).astype(int)]
    # contribution from higher bin
    froll = frac * irf[np.fmod(fwd + n, n).astype(int)]
    return broll + froll

def exp_decay(time, tau):
    return np.exp(-time / tau)

def nnls_conv_irf(x_in, irf, params, y_in):
    irf_shift, *tau0 = params
    shifted_irf = shift_irf(irf, irf_shift)
    decays = [convolve(shifted_irf, exp_decay(x_in, tau)) for tau in tau0]
    decays.append(np.ones_like(x_in))  # Adding a constant offset term
    A = np.vstack(decays).T
    x_out, _ = nnls(A, y_in)
    return A, x_out, np.dot(A, x_out)

def residual_function(params, x_in, irf, y_in):
    _, _, y_fit = nnls_conv_irf(x_in, irf, params, y_in)
    residuals = (y_fit - y_in) / np.sqrt(y_fit + 1)
    chi2 = np.sum(residuals ** 2)
    chi2_reduced = chi2 / (len(x_in))
    return chi2_reduced

def reconvolution_fit(data, exp_num=1, tau_bounds=None, maxiter=1000,
                      disp=False, workers=1):
    '''
    perform reconvolution fit of TCSPC data using differential evolution.
    note that this is a stochastic method; sometimes (especially when
    doing monoexponential fits, in my experience) the output might be
    nonsensical. for now just try running it again.
    '''
    x     = data["Time(s)"].to_numpy()
    y     = data["Emitted"].to_numpy()
    irf   = data["IRF"].to_numpy()
    if tau_bounds is None:
        tau_bounds = [(1e-12, 1.)] * exp_num
    tau_bounds = [(max(bound[0], 1e-12),
                   max(bound[1], 1e-12)) for bound in tau_bounds]
    param_bounds = [(-len(x)/2., len(x)/2.)]
    param_bounds.extend(tau_bounds)

    # Global minimization using differential evolution
    result = differential_evolution(residual_function, bounds=param_bounds, 
            args=(x, irf, y), strategy='best1bin', maxiter=maxiter,
            tol=1e-12, popsize = 5*exp_num*3, polish = True,
            workers = workers, disp = disp)
    best_popt = result.x
    irf_shift_opt = best_popt[0]
    tau_opt = best_popt[1:]
    # pull out the best-fit amplitudes from the NNLS
    A, x_out, y_fit = nnls_conv_irf(x, irf, best_popt, y)
    amplitudes = x_out[:-1]
    offset = x_out[-1]
    chi2_reduced = residual_function(best_popt, x, irf, y)
    # Sort tau and amplitudes based on tau
    tau_opt, amplitudes = zip(*sorted(zip(tau_opt, amplitudes),
                                      key=lambda x: x[0]))
    return tau_opt, amplitudes, irf_shift_opt, offset, chi2_reduced

def multi_fit(filename, nmax):
    df = get_histogram(filename)
    x = df['Time(s)']
    y = df['Emitted']
    # make a dict of the fits
    od = {}
    od['Time(s)'] = x
    od['Emitted'] = y
    logfile = os.path.splitext(filename)[0] + "_fit_log.txt"
    with open(logfile, "w") as f:
        for n_exp in range(1, nmax + 1):
            print(f"Fitting histogram from {filename} "
            f"with {n_exp} exponentials.")
            tau_opt, amps_opt, shift_opt, offset, rchi2 = reconvolution_fit(
                df,exp_num=n_exp)
            f.write(f"Fitting with n_exp = {n_exp}:\n")
            f.write(f"tau = {tau_opt}\n")
            f.write(f"amps = {amps_opt}\n")
            f.write(f"irf_shift = {shift_opt}\n")
            f.write(f"offset = {offset}\n")
            f.write(f"reduced chi^2 = {rchi2}\n") 
            f.write("\n")
            title = [f"n_exp = {n_exp}"]
            taustr = (', ').join([sci_format(t) for t in tau_opt])
            ampstr = (', ').join([f"{a:.4f}" for a in amps_opt])
            title.append(f"taus = {taustr}")
            title.append(f"amps = {ampstr}")
            title.append(f"irf_shift = {shift_opt:.2}")
            title.append(f"offset = {offset:.2}")
            title.append(f"reduced_chi^2 = {rchi2:.2}")
            tstr = ('\n').join(title)
            opt = np.zeros_like(x)
            for j in range(n_exp):
                opt += amps_opt[j] * exp_decay(x, tau_opt[j])
            od[f"Decays_{n_exp}_exp"] = opt
            shifted_irf = shift_irf(df['IRF'], shift_opt)
            od[f"Shifted_IRF_{n_exp}_exp"] = shifted_irf
            y_fit = convolve(opt, shifted_irf)
            od[f"Fit_{n_exp}_exp"] = y_fit
            ymax = np.max(df['Emitted'])
            residuals = (y_fit - y) / np.sqrt(y_fit + 1)
            fig, ax = plt.subplots(2, 1,
                        gridspec_kw={'height_ratios': [3, 1]},
                        layout="constrained",
                        sharex=True, figsize=(12, 10))
            ax[0] = fix_y(ax[0], ymax)
            ax[0].plot(x, df['Emitted'], label="Emitted", alpha=0.5, lw=2.)
            ax[0].plot(x, shifted_irf * (ymax / np.max(shifted_irf)),
                       label="IRF (scaled + shifted)", alpha=0.5, lw=2.)
            ax[0].plot(x, y_fit,
                       label=r'Fit ($ n_{\text{exp}} ' + f" = {n_exp} $)", lw=3.)
            ax[0].legend()
            ax[1] = fix_x(x, 4, ax[1])
            ax[1].plot(x, residuals, label='Residuals', alpha=0.5, lw=2.)
            ax[1].set_ylabel('Residuals')
            ax[1].axhline(y=0, color='grey', linestyle='--')
            ax[1].legend()
            for axis in ax:
                axis.grid(visible=True, which='major', axis='both',
                          color='k', alpha=0.25, linestyle='--', lw=0.5)
            outfile = os.path.splitext(filename)[0] + f"n_exp_{n_exp}_ml_fit.pdf"
            fig.suptitle(tstr)
            fig.savefig(outfile)
            plt.close()
    fdf = pd.DataFrame(od)
    fdf_file = os.path.splitext(filename)[0] + "_fits.csv"
    fdf.to_csv(fdf_file)
