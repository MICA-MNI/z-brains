import numpy as np
def compute_fmri_features(data, tr, f_low=0.01, f_high=0.08):
    """
    compute multiple fMRI features: rmssd, timescales, alff, falff.

    parameters
    ----------
    data : array-like, shape (n_timepoints, n_vertices)
        time × vertex (or region) matrix.
        example: data[t, v] is the signal at time t for vertex v.
    tr : float
        repetition time in seconds.
    f_low : float, optional
        lower bound of frequency band for alff/falff (hz), default 0.01.
    f_high : float, optional
        upper bound of frequency band for alff/falff (hz), default 0.08.
    returns
    -------
    """
    rmssd = mica_func_rmssd(data)
    acf, timescales = mica_func_timescales(data, tr)
    alff, falff = mica_func_alff(data, tr, f_low, f_high)
    return {
        "rmssd": rmssd,
        "timescales": timescales,
        "alff": alff,
        "falff": falff
    }


def mica_func_rmssd(data):
    """
    compute rmssd (root mean square of successive differences).

    parameters
    ----------
    data : array-like, shape (n_timepoints, n_vertices)
        time × vertex (or region) matrix.
        example: data[t, v] is the signal at time t for vertex v.

    returns
    -------
    rmssd : ndarray, shape (n_vertices,)
        rmssd value for each vertex.

    usage
    -----
    >>> import numpy as np
    >>> rmssd = mica_func_rmssd(data)
    >>> np.savetxt("rmssd.csv",  alff.reshape(1, -1), delimiter=",")
    """
    data = np.asarray(data, dtype=float)
    n_timepoints, n_vertices = data.shape

    rmssd = np.zeros(n_vertices, dtype=float)

    # mean across time for each vertex (axis=0)
    mean_ts = data.mean(axis=0)  # shape: (n_vertices,)

    for j in range(n_vertices):
        centered = data[:, j] - mean_ts[j]
        differences_squared = np.diff(centered) ** 2
        rmssd[j] = np.sqrt(np.sum(differences_squared)) / (n_timepoints - 1)

    return rmssd

def mica_func_timescales(data, tr):
    """
    compute intrinsic timescales using the autocorrelation function.

    parameters
    ----------
    data : array-like, shape (n_timepoints, n_vertices)
        time × vertex (or region) matrix.
        example: data[t, v] is the signal at time t for vertex v.
    tr : float
        repetition time in seconds.

    returns
    -------
    acf : ndarray, shape (n_timepoints-1, n_vertices)
        autocorrelation function at positive lags (excluding lag 0).
    timescales : ndarray, shape (n_vertices,)
        intrinsic timescale for each vertex, in seconds.

    usage
    -----
    >>> import numpy as np
    >>> acf, timescales = mica_func_timescales(data, tr)
    >>> np.savetxt("timescales.csv", timescales.reshape(1, -1), delimiter=",")
    """
    data = np.asarray(data, dtype=float)
    n_timepoints, n_vertex = data.shape

    acf = np.zeros((n_timepoints - 1, n_vertex), dtype=float)
    timescales = np.zeros(n_vertex, dtype=float)

    for i in range(n_vertex):
        signal = data[:, i].astype(float)

        # full autocorrelation
        acf_full = np.correlate(signal, signal, mode='full')

        # normalize so that lag 0 = 1
        zero_lag = acf_full[n_timepoints - 1]
        if zero_lag != 0:
            acf_full = acf_full / zero_lag

        # positive lags only (exclude lag 0)
        acf_i = acf_full[n_timepoints:]  # length = n_timepoints - 1
        acf[:, i] = acf_i

        # sum until the first negative value
        for t in range(n_timepoints - 1):
            if acf_i[t] >= 0:
                timescales[i] += acf_i[t]
            else:
                break

    timescales = timescales * tr
    return acf, timescales

def mica_func_alff(data, tr, f_low, f_high):
    """
    compute alff and falff within a custom frequency band.

    parameters
    ----------
    data : array-like, shape (n_timepoints, n_vertices)
        time × vertex (or region) matrix.
        example: data[t, v] is the signal at time t for vertex v.
    tr : float
        repetition time in seconds.
    f_low : float
        lower bound of frequency band (hz), e.g., 0.01.
    f_high : float
        upper bound of frequency band (hz), e.g., 0.08.

    returns
    -------
    alff : ndarray, shape (n_vertices,)
        amplitude of low-frequency fluctuations for each vertex.
    falff : ndarray, shape (n_vertices,)
        fractional alff for each vertex (alff / total amplitude).

    usage
    -----
    >>> import numpy as np
    >>> alff, falff = mica_func_alff(data, tr, f_low=0.01, f_high=0.08)
    >>> np.savetxt("alff.csv",  alff.reshape(1, -1), delimiter=",")
    >>> np.savetxt("falff.csv", falff.reshape(1, -1), delimiter=",")
    """
    data = np.asarray(data, dtype=float)
    n_timepoints, n_vertices = data.shape

    # detrend each vertex (remove best-fitting straight line)
    x = np.arange(n_timepoints)
    data_detrend = np.zeros_like(data)
    for i in range(n_vertices):
        y = data[:, i]
        a, b = np.polyfit(x, y, 1)   # linear fit y ≈ a*x + b
        trend = a * x + b
        data_detrend[:, i] = y - trend

    # one-sided fft along time dimension (non-negative frequencies only)
    y_fft = np.fft.rfft(data_detrend, axis=0)
    freqs = np.fft.rfftfreq(n_timepoints, d=tr)
    amp = np.abs(y_fft)

    # mask for custom low-frequency band
    mask_low = (freqs >= f_low) & (freqs <= f_high)

    # compute alff and falff
    alff = amp[mask_low, :].sum(axis=0)                # shape: (n_vertices,)
    total_amp = amp.sum(axis=0)                        # shape: (n_vertices,)

    falff = np.zeros_like(alff)
    nonzero = total_amp > 0
    falff[nonzero] = alff[nonzero] / total_amp[nonzero]
    falff[~nonzero] = np.nan

    return alff, falff