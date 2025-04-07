import os
import sympy
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def get_histogram(filename):
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
    emissive_str = lines[1].split(" ")
    emissive = [True if s == "T" else False for s in emissive_str]
    str_array = [line.split(" ") for line in lines[2:]]
    bins = np.array([row[0] for row in str_array]).astype(float)
    counts = np.array([row[1:] for row in str_array]).astype(float)
    emissive_counts = counts[:, emissive[1:]]
    if emissive_counts.shape[1] > 1:
        sum_emissive = np.sum(emissive_counts, axis=0)
    else:
        sum_emissive = emissive_counts
    return labels, bins, counts, sum_emissive

def plot_all(labels, bins, counts, outfile):
    '''
    big plot with all binned decay pathways plotted, including
    ones which would be invisible to a real detector.
    '''
    fig, ax = plt.subplots(figsize=(12,8))
    for i, l in enumerate(labels):
        if np.sum(counts[:, i]) > 0:
            plt.plot(bins, counts[:, i], label=l)
    plt.grid(visible=True)
    ax.set_xlabel("time")
    ax.set_ylabel("counts (unnormalised)")
    ax.legend()
    ax.set_yscale('log')
    fig.tight_layout()
    plt.savefig(outfile)
    plt.close()

'''
below is all taken from my old aggregate code
'''

def Convol(x, h):
    X = np.fft.fft(x)
    H = np.fft.fft(h)
    return np.real(np.fft.ifft(X * H))

def exp_model(t, *args):
    '''
    given a set of amplitudes and time constants as
    (a_1, ..., a_n, tau_1, ..., tau_n), return the
    multiexponential curve as a function of t
    '''
    if len(args) == 0 or len(args) % 2 != 0:
        raise ValueError("exp_model: incorrect number of arguments")
    y = np.zeros(t.size)
    n_exp = (len(args)) // 2
    for i in range(n_exp):
        y += args[i] * np.exp(-t / float(args[i + (n_exp)]))
    return y

def reconv(X, *args):
    '''
    return the multiexponential reconvolved with the IRF.
    this is ugly internally as described lower in do_fit;
    basically because of some limitation in scipy.curve_fit
    X should be a tuple of t, the IRF, and the time constants.
    args should be a_i, irf_shift
    '''
    t, irf, taus = X
    ymodel = np.zeros(t.size)
    irf_interp = np.interp(t, t - args[-1], irf)
    irf_reshaped_norm = irf_interp / np.sum(irf_interp)
    for i in range(len(args) - 1):
        ymodel += (args[i] * np.exp(-(t) / taus[i]))
    z=Convol(ymodel,irf_reshaped_norm)
    return z

def lifetimes(n, names, bv, err, covar):
    '''
    given a set of n exponentials, construct the expressions for
    tau_amp and tau_int, calculate them with the best values `bv`,
    calculate the error on those, and return a dict with all the
    relevant information.
    '''
    strings = [[], [], []]
    # define ai, taui for i = 1, n
    sympy.symbols('tau:{:d}, a:{:d}'.format(n, n))
    # build up the expressions for ai, ai * taui, ai * taui^2
    for i in range(1, n + 1):
        for j in range(3):
            strings[j].append('a{:d} * tau{:d}**{:d}'.format(i, i, j))
    # turn the lists of strings into the relevant sympy expressions
    joined = [' + '.join(s) for s in strings]
    tau = [sympy.sympify(j, evaluate=False) for j in joined]
    # we need UnevaluatedExpr here otherwise sympy cancels the a1 for
    # the monoexponential fit and we never get out its value or error
    tau_amp = sympy.UnevaluatedExpr(tau[1]) / sympy.UnevaluatedExpr(tau[0])
    tau_int = sympy.UnevaluatedExpr(tau[2]) / sympy.UnevaluatedExpr(tau[1])
    # now start on relating these expressions to the fitted parameters
    j_amp = np.zeros(len(names))
    j_int = np.zeros(len(names))
    var = list(tau_int.free_symbols) # same for both amp and int
    tau_amp = tau_amp.doit()
    tau_int = tau_int.doit()
    # generate a list of tuples which tell sympy the values to substitute in
    repl = [(var[i], bv[str(var[i])]) for i in range(len(var))]
    # we're gonna return a dict which we turn into a pandas dataframe
    # then compare to find how many exponents gave the best fit
    d = {'n_exp': n}
    for i in range(len(var)):
        # dict key and index comparison require string representation
        s = str(var[i])
        # build up the dict as we go
        d[s] = bv[s]
        d[s + '_err'] = err[s]
        '''
        sympy doesn't order the variables ai, taui, but returns them as a set.
        they are ordered in curve_fit though - whatever order we put them in in,
        the covariance matrix etc is ordered the same way. so use `names` to find
        the right index to put the derivative in and use that.
        note that this also leaves the indices corresponding to x0 and y0 = 0,
        wherever they are in the list, so we don't need to worry about them.
        '''
        ind = np.nonzero(np.array(names) == s)[0][0]
        j_amp[ind] = sympy.diff(tau_amp, var[i]).subs(repl)
        j_int[ind] = sympy.diff(tau_int, var[i]).subs(repl)
    m_amp = np.matmul(j_amp, covar)
    m_int = np.matmul(j_int, covar)
    tau_amp_err = np.sqrt(np.matmul(m_amp, j_amp.transpose()))
    tau_int_err = np.sqrt(np.matmul(m_int, j_int.transpose()))
    d['tau_amp'] = tau_amp.subs(repl)
    d['tau_amp_err'] = tau_amp_err
    d['tau_int'] = tau_int.subs(repl)
    d['tau_int_err'] = tau_int_err
    print("n_exp = {:d}".format(n))
    print("tau_amp = {} +/- {} ps".format(tau_amp.subs(repl), tau_amp_err))
    print("tau_int = {} +/- {} ps".format(tau_int.subs(repl), tau_int_err))
    return d

def do_fit(filename, tau_init, irf_file=None,
        exp=False, pw=0.25, pm=0.4, time_unit="ns"):
    path = os.path.dirname(filename)
    (fluence, ext) = os.path.splitext(os.path.basename(filename))
    
    '''
    import the decay data. 
    if it's experimental it should be a two-column file with the first column
    being the time and the second either being counts or "normalised" counts.
    otherwise it's the histogram output by the fortran code: time, ann, pool, pq, q
    '''
    labels, bins, all_counts, ec = get_histogram(filename)
    ecn = ec / np.max(ec)
    xyn = np.column_stack((bins, ec, ecn))
    max_count_time = xyn[np.argmax(ec), 0]
    max_time = np.max(bins)
    
    cutoff = max_count_time

    # errors for each count
    if np.max(xyn[:, 1]) > 1.:
        max_count = np.max(xyn[:, 1])
    else:
        print("Warning - assuming max_count = 10000. bin errors might be wrong")
        max_count = 10000. # arbitrary!
    sigma = np.zeros(xyn[:, 1].size)
    for i in range(len(xyn[:, 1])):
        if (xyn[i, 1] == 0.):
            count = 1.
        else:
            count = xyn[i, 1]
        sigma[i] = np.sqrt(1. / count + 1. / max_count)
        
    n_exp = len(tau_init)
    p0 = [*[1./n_exp for _ in range(n_exp)], *tau_init]
    names = [] # these will be needed to construct a dataframe later
    for i in range(n_exp):
        names.append("a{:d}".format(i + 1))
    for i in range(n_exp):
        names.append("tau{:d}".format(i + 1))
    # bounds for each of the time constants
    lbs = tuple([0. for _ in range(len(p0))])
    ubs = tuple([np.inf for _ in range(len(p0))])
    bounds = [lbs, ubs]

    tail = xyn[xyn[:, 0] >= cutoff]
    tail_sigma = sigma[xyn[:, 0] >= cutoff]
    x = tail[:, 0] - cutoff
    y = tail[:, 2]
    tail_popt, tail_pcov = curve_fit(exp_model, x, y, p0=p0,
            sigma=tail_sigma, bounds=bounds)
    print("popt tail: ",tail_popt)
    print(tail_pcov)
    tail_err = np.sqrt(np.diag(tail_pcov))
    print("errors:", tail_err)
    
    best_t = list(tail_popt[len(tail_popt)//2:])
    print("Time constant(s) from tail fit = ", best_t)
    best_a = list(tail_popt[:len(tail_popt)//2])
    print("Amplitude(s) from tail fit = ", best_a)
    # need this to plot the tail fit later
    bf = exp_model(xyn[:, 0], *tail_popt)
    
    '''
    now do the same for the IRF
    '''
    irf_norm = np.zeros(ec.size)
    # NB: update this. the fortran should output a pulse file
    # and then we read that in
    sig = pw / 2.355
    irf_gen = ((1 / (sig * np.sqrt(2. * np.pi))) *
            np.exp(-(xyn[:, 0] - pm)**2 / (np.sqrt(2.) * sig)**2))
    irf_norm = irf_gen / np.max(irf_gen)
    
    # fit tail with IRF
    fig, ax = plt.subplots(figsize=(12,8))
    plt.plot(x, y, ls='--', marker='o', label='Decays')
    plt.plot(x, exp_model(x, *tail_popt),
            label=r'Fit: $ \tau_i = $' + f"{best_t}")
    plt.plot(xyn[:, 0] - cutoff, irf_norm, label='IRF')
    plt.legend(fontsize=32)
    plt.grid()
    ax.set_yscale('log')
    ax.set_ylim([1e-5, 1.1])
    ax.set_ylabel("Counts (normalised)")
    ax.set_xlabel("time")
    plt.tight_layout()
    plt.savefig(os.path.join(path, f"{fluence}_tail_fit_{n_exp}.pdf"))
    plt.close()

    """
    now we need to do something horrible!
    generate an empty array with the same length as x and irf, and fill
    the first n_exp elements with our time constants. 
    we do this because we want to keep them fixed in the subsequent fit, but:
      - you can't just pass (X, irf, *best_t) because
        the arrays have to have the same shape (?)
      - you can't wrap it in lambda X, *best_t because
        then the reconvolution function is passed
        with two extra arguments, for some reason
    This second point is something internal to curve_fit;
    just doing f = lambda X, *best_t: 
    print(len((X, *best_t, *best_a, x0, irf_shift))) returns 5, but when you
    do that in curve_fit it returns 7. 
    no idea why. bug? something to do with self?
    """
    
    taus = np.zeros(len(xyn[:, 0]))
    for i in range(n_exp):
        taus[i] = best_t[i]

    irf_shift = 0.0
    X = (xyn[:, 0], irf_norm, taus)
    p0 = [*best_a, irf_shift]
    lbs = tuple([0. for _ in range(len(best_a))] + [-max_time])
    ubs = tuple([np.inf for _ in range(len(best_a))] + [max_time])
    bounds = [lbs, ubs]
    popt, pcov = curve_fit(reconv,
            X, xyn[:, 2], p0=p0, sigma=sigma, bounds=bounds)
    bf = reconv(X, *popt)

    print("best fit for amps, irf_shift: ", popt)
    print("cov: ", pcov)

    # now set up the lifetime calculations
    err = np.sqrt(np.diag(pcov))
    print("err: ", err)
    best_values = dict(zip(names, np.concatenate((popt[:n_exp], best_t))))
    errors = dict(zip(names, np.concatenate((err[:n_exp], tail_err[n_exp:]))))
    # amplitudes and time constants are fitted separately
    # so their covariance is necessarily zero
    cov = np.block([
        [pcov[:n_exp, :n_exp], np.zeros((n_exp, n_exp))],
        [np.zeros((n_exp, n_exp)), tail_pcov[n_exp:, n_exp:]]
    ])

    d = lifetimes(n_exp, names, best_values, errors, cov)
    d["n_exp"] = n_exp
    d["irf_shift"] = popt[-1]
    d["irf_shift_err"] = err[-1]
    d["cutoff"] = cutoff
    print(d)

    fig, axs = plt.subplots(2, 1, figsize=(12,8))
    plt.suptitle("{}: ".format(fluence) + r'$ \tau_{\text{amp.}} = $'
            + "{:4.2f} +/- {:4.2f} ns".format(d["tau_amp"], d["tau_amp_err"]))
    axs[0].plot(xyn[:, 0], xyn[:, 2], ls='--', marker='o', label='Decays')
    axs[0].plot(xyn[:, 0], bf, label='fit')
    plot_file = os.path.join(path, f"{fluence}_reconv_{n_exp}.pdf")
    axs[0].legend()
    axs[0].grid(True)
    axs[1].grid(True)
    axs[0].set_ylim([0., 1.1])
    axs[0].set_ylabel("Counts")
    axs[0].set_xlabel("Time (ns)")
    axs[0].set_xlim([-1., 10.])
    axs[1].plot(xyn[:, 0], xyn[:, 2], ls='--', marker='o', label='Decays')
    axs[1].plot(xyn[:, 0], reconv(X, *popt), label='fit')
    axs[1].set_yscale('log')
    axs[1].set_ylim([1e-2, 1.5])
    axs[1].set_ylabel("Counts")
    axs[1].set_xlabel("Time (ns)")
    axs[1].set_xlim([-1., 10.])
    fig.tight_layout()
    plt.savefig(plot_file)
    plt.close()
    df = pd.DataFrame(d, index=[0])
    df_file = os.path.join(f"{path}", f"{fluence}_fit_{n_exp}.csv")
    df.to_csv(df_file)

    count_file = os.path.join(f"{path}", f"{fluence}_norm_counts.txt")
    np.savetxt(count_file, np.column_stack((xyn[:, 0], xyn[:, 2])))
    fit_file = os.path.join(f"{path}", f"{fluence}_fit_xy_{n_exp}.txt")
    np.savetxt(fit_file, np.column_stack((xyn[:, 0], bf)))

    return (d, np.column_stack((xyn[:, 0], xyn[:, 2], bf)))
