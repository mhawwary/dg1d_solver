from numpy.fft import rfft, rfftfreq, fft, fftfreq
from numpy import size, shape, mean, conj, hanning

#-Compute FFT for real input data
# takes advantage of having negative frequency components
# as the conjugate of the positive frequency componenets
def compute_rfft(u_data_, sample_spacing):
    sample_size = size(u_data_)
    # fftsize: the size of half the spectra with positive frequencies
    if sample_size%2==0:  # even
        fftsize = int((sample_size/2)+1)
    elif sample_size%2==1: # odd
        fftsize  = int((sample_size+1)/2)
    u_fft = rfft(u_data_)/fftsize     # fft
    freq  = rfftfreq(sample_size,sample_spacing) # frequencies in cycles per sec
    return freq, u_fft

#-Compute FFT for complex input data
def compute_cfft(u_data_, sample_spacing):
    sample_size = size(u_data_);
    # fftsize+1: the size of half the spectra with positive frequencies
    if sample_size%2==0:
        fftsize = int(sample_size/2)
    elif sample_size%2==1:
        fftsize = int((sample_size-1)/2)
    u_fft_ = fft(u_data_)/fftsize    # fft
    freq_  = fftfreq(sample_size,sample_spacing) # frequencies in cycles per sec
    # taking half the spectra + u_fft_[n/2] where the Nyquist freq exist
    u_fft = u_fft_[0:fftsize+1]   # fft_size = Nt/2 for Nt even, so the size of this array is fftsize+1
    freq  = freq_[0:fftsize+1]
    freq[-1]=freq[-1]*(-1)**(sample_size+1); # adjusting the sign of this frequency

    return freq, u_fft

# This function computes the ensemble averaged fft
# for a set of data for KE spectra
# it outputs freq as the integer wavenumbers
def compute_ensemble_fft(u_data_,sample_spacing):
    sample_size = size(u_data_[0,:])   # for one fft data array
    sample_spacing = 1./sample_size    # for integer wavenumbers
    aver_size = size(u_data_[:,0])     # n of fft data arrays
    if sample_size%2==0: # even
        fftsize = int((sample_size/2)+1) # the size of half the spectra with positive frequencies
    elif sample_size%2==1: # odd
        fftsize = int((sample_size+1)/2)
    freq, u_fft = compute_rfft(u_data_[0,:],sample_spacing) # freq will be integer wavenumbers
    u_amp = abs(u_fft)  # amplitude
    for i in range(1,aver_size):
        u_fft = compute_rfft(u_data_[i,:],sample_spacing)
        u_amp += abs(u_fft)
    u_amp = u_amp / (aver_size*fftsize)
    KEnerg = 0.5*u_amp**2.0

    return freq, u_amp, KEnerg

# Power Spectral Density for real inputs
# needs to include shifting and averaging
def compute_psd(u_data_, sample_spacing):
    sample_size = size(u_data_)
    if sample_size%2==0:  # even
        fftsize = int((sample_size/2)+1) # size of half the spectra
    elif sample_size%2==1: # odd
        fftsize = int((sample_size+1)/2)
    u_data_ = u_data_ - mean(u_data_)
    freq = fft.rfftfreq(sample_size, sample_spacing)
    u_fft = fft.rfft(u_data_)/fftsize
    u_psd = u_fft * conj(u_fft)

    return freq, u_psd

