import numpy as np
import sys
import wavelet as wavelet
from helpers import find, ar1_spectrum


########################################################
def cwt_results(data, dt, lag1, dj=1/12, s0=-1, J=-1, mother='morlet', slevel=0.95):
    """
    Main function to get wavelet spectrum
    
    Parameters
    ----------
    data : numpy.ndarray
        Input signal array with length=N
    dt : float
        Sampling interval.
    dj : float, optional
        Spacing between discrete scales. Default value is 1/12.
        Smaller values will result in better scale resolution, but
        slower calculation and plot.
    s0 : float, optional
        Smallest scale of the wavelet. Default value is 2*dt.
    J : float, optional
        Number of scales less one. Scales range from s0 up to
        s0 * 2**(J * dj), which gives a total of (J + 1) scales.
        Default is J = (log2(N * dt / so)) / dj.
    mother : Wavelet class, or string
        Mother wavelet class. Default is Morlet wavelet.
    freqs : numpy.ndarray, optional
        Custom frequencies to use instead of the ones corresponding
        to the scales described above. Corresponding scales are
        calculated using the wavelet Fourier wavelength.

    Returns
    -------
    freq : array like
        Vector of Fourier frequencies (1 / time units) that
        corresponds to the wavelet scales.
    period : array like
        Vector of Fourier periods (1 / freq).
    scale : numpy.ndarray
        Vector of scale indices given by sj = s0 * 2**(j * dj),
        j={0, 1, ..., J}.  
    wave : numpy.ndarray
        Wavelet transform according to the selected mother wavelet.
        Has (J+1) x N dimensions.
    iwave : numpy.ndarray
        Inverse wavelet transform according to the selected mother wavelet.
        Should be as close as possible to the input data.
    power : numpy.ndarray
        The power spectrum = wave**2
    powerc : numpy.ndarray
        The power spectrum with bias correction (Liu et al., 2007).
    signif : array like
        Significance levels as a function of scale or period (and scaled with variance of signal, =1 for standardized data).
        A regular chi-square test is performed according to Torrence and Compo (1998) equation 18.
    coi : numpy.ndarray
        Returns the cone of influence, which is a vector of N points containing the maximum Fourier period of useful
        information at that particular time. Periods greater than those are subject to edge effects.
    var_exp_wave : float
        variance of the wavelet (Torrence and Compo, 1998, eq. (14)) / variance of the data * 100. 
    rmse_wave : float
        root mean square error between the inverse wavelet and the data 

    """
    #info from input data-------
    N=len(data)
    var_data=np.var(data,ddof=0) #variance of input data. =1 when input data is standardized

    #wavelet transform--------
    #get the mother class (default 'morlet')
    mother = wavelet.check_parameter_wavelet(mother)
    #wavelet transform (mandatory remove mean but not necessary a standarization)
    #wave, scale, freq, coi, fft, fft_freq = wavelet.cwt(data, dt, dj, s0, J, mother)
    wave, scale, freq, coi = wavelet.cwt(data, dt, dj, s0, J, mother)
    #                                                         
    #the power spectrum 
    power = (np.abs(wave))**2 #dimensionless for standardized data, or units of variance for non-standardized data
    #
    #the Fourier equivalent periods for each wavelet scale
    period = 1./freq
    # 
    #The power is significant when the ratio power / signif > 1
    #signif is scaled with the variance (if data is standardized, signif=1 and dimensionless)
    signif,_ = wavelet.significance(var_data, dt, scale, 0, lag1, significance_level=slevel, wavelet=mother)
    signif = np.ones([1, N]) * signif[:, None]; signif = power / signif

    #bias correction--------
    #- Original bias correction = divide the power with the scale (Liu et al., 2007)
    #- Our proposal for bias correction = divide the power with scale and cdelta (cte) and multiply by dt
    #                 so the power spectrum has units of variance of data, or dimensionless if data is standardized.
    #                 It represents the local contributions of the power to the variance and is independent of the resolution dj.
    #                 So this is a standardized correction more consistent with energy (variance of original or standardized data)
    #- High freq peaks are improved, but the same peaks are significant in the cases with and without bias corrections.
    #- The significance ratio "signif" is the same when using or not bias correction
    #- It highlights better the dominant forcing frequency.
    powerc=(power.T/scale).T*dt/mother.cdelta #powerc is dimensionless (for standardized data) or has units of variance (for non-standardized data)
    #signifc = signif/scale*dt/mother.cdelta
    #signif = power/signif = powerc/signifc

    #accuracy of wavelet transform---------
    #check the if wavelet has the 100% of original variance
    #a good wavelet set of parameters should capture almost 100% of the data variance
    #variance of wavelet spectrum: as in Torrence and Compo (1998) eq. (14)
    var_wave=dj/N*np.sum(powerc) #ideally should be = 1 if data is standardized, or = variance of data for non-standardized data
    var_exp_wave=np.round(var_wave/var_data*100,1) #percentage of original variance retain in the wavelet
    #
    #check the root mean square error between original data and inverse wavelet
    iwave = wavelet.icwt(wave, scale, dt, dj, mother) #inverse transform (units of data or dimensionless)
    rmse_wave=(np.sum((iwave-data)**2)/N)**.5 #units of data or dimensionless

    #save results-----------
    #results from wavelet
    result_wavelet = {'freq':freq,'period':period,'scale':scale,
                      'wave':wave,'iwave':iwave,'power':power,'powerc':powerc,'signif':signif,'coi':coi,
                      'var_exp_wave':var_exp_wave,'rmse_wave':rmse_wave} 
    #
    return result_wavelet
########################################################


########################################################
def cwt_filter(period, scale, wave, dt, dj, mother, filter_type, cutoff_periods):
    """
    Perfomr a wavelet filter = inverse wavelet transform for some scales or periods
    
    Parameters
    ----------
    period : array like
        Vector of Fourier periods (1 / freq) equivalent to scales.
    scale : numpy.ndarray
        Vector of scale indices given by sj = s0 * 2**(j * dj),
        j={0, 1, ..., J}.   
    wave : numpy.ndarray
        Wavelet transform according to the selected mother wavelet.
        Has (J+1) x N dimensions.
    dt : float
        Sampling interval.
    dj : float, optional
        Spacing between discrete scales. Default value is 1/12.
        Smaller values will result in better scale resolution, but
        slower calculation and plot.
    mother : Wavelet class, or string
        Mother wavelet class. Default is Morlet wavelet.
    filter_type : str
                   Type of wavelet filter to perform. 
                   "low-pass" or "high-pass" or "band-pass".
    cutoff_periods : list
                      List used to draw cut-off periods of data filter with wavelets.
                      [ref_period] = for low-pass or high-pass filters.
                      [period1, period2] = for band-pass filter.

    Returns
    -------
    iwave_range : array like
        Filtered time series according to the filter_type selected.
    """
    #
    filter_types=["high-pass","low-pass","band-pass"]
    if filter_type in filter_types:
        if filter_type=="high-pass":
            sel = find(period <= cutoff_periods[0]) #find periods < cutoff_periods[0]    
        elif filter_type=="low-pass":
            sel = find(period >= cutoff_periods[0]) #find periods > cutoff_periods[0]                       
        else:
            sel = find((period >= cutoff_periods[0]) & (period <= cutoff_periods[1])) #find periods in range of interest
    else:
        sys.exit(f"use a filter_type name from this list: {filter_types}")
    #get the mother class (default 'morlet')
    mother = wavelet.check_parameter_wavelet(mother)
    #As of Torrence and Compo (1998) eq. (29)
    cte=np.abs(dj*np.sqrt(dt)/(mother.cdelta*mother.psi(0)))
    iwave_range=(np.real(wave).T/np.sqrt(scale)).T
    iwave_range=cte*iwave_range[sel, :].sum(axis=0) #reconstruct time serie using the above periods
    #
    return iwave_range
########################################################