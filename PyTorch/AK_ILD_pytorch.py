import torch
import numpy as np
from scipy.fft import rfft,irfft
from scipy.signal import lfilter


def clc_e_ILD(P, f_c, fs, f_lim, Nmax, n, C, nfft,ILD_ref):
    P = torch.permute(P, (0, 2, 1))
    x_l = torch.squeeze(P[:, 0, :]) # space x left/right x time (361, 2, 768)
    x_l = x_l.transpose(0, 1)
    x_r = torch.squeeze(P[:, 1, :])
    x_r = x_r.transpose(0, 1)
    # get length and channels
    N = x_l.shape[0]  # ir length
    c = x_l.shape[1]  # ir num. of channels
    

    # Filter x with filter bank and calculate error in each band
    ILD = torch.zeros((n, c))
    C_db_tmp = 10 * torch.log10(torch.abs(C))
    C_db = torch.zeros_like(C_db_tmp)
    # ------------------------------------- 2. filter signals and get ERB error
    for k in range(n):
        # find point where filter decayed for 60 dB
        C_db = C_db_tmp[:, k]
        #Nmin = torch.argmax(torch.flip(C_db, dims=[0]) > torch.max(C_db) - 60) + 1
        #Nmin = 2 ** torch.ceil(torch.log2(Nmax - Nmin + 1))
        # needed filter length at current center frequency
        #Ncur = 2 ** torch.ceil(torch.log2(1 / f_c[k] * 4 * fs))
        #Ncur = int(max(Ncur, Nmin, N))  # consider replacing with Ncur = N?
        
        Ncur = N

        # get indices of upper and lower frequency limit
        f_lim_id = [round(f_lim[0] / (fs / Ncur)) + 1, round(f_lim[1] / (fs / Ncur)) + 1]
        f_lim_id[0] = max(f_lim_id[0], 2)  # make sure we don't use the 0 Hz bin

        # filter input signals
        t11 = torch.abs(torch.fft.fft(C[:, k], Ncur, dim=0)).reshape(-1, 1)
        t1 = t11.repeat(1, c)
        t2 = torch.abs(torch.fft.fft(x_l, Ncur, dim=0)) ** 2
        X_L = t1 * t2

        t2 = torch.abs(torch.fft.fft(x_r, Ncur, dim=0)) ** 2
        X_R = t1 * t2

        # get energy
        X_L = torch.sum(X_L[f_lim_id[0]:f_lim_id[1], :], dim=0)
        X_R = torch.sum(X_R[f_lim_id[0]:f_lim_id[1], :], dim=0)

        # get ERB error
        ILD[k, :] = 10 * torch.log10(torch.abs(X_L) / torch.abs(X_R))
        # ILD[k, :] = torch.abs(X_L) / torch.abs(X_R)
        # ILD[k, :] = X_L / X_R
    
    e_ILD_omega    = torch.mean(torch.abs(ILD_ref - ILD),dim=0) #avarage over frequencies bands
    return e_ILD_omega


def ERBFilterBank(x, fcoefs=None):
    if x.ndim < 1:
        raise ValueError("Syntax: output_array = ERBFilterBank(input_vector[, fcoefs]);")

    if fcoefs is None:
        fcoefs = MakeERBFilters(22050, 64, 100)

    if fcoefs.shape[1] != 10:
        raise ValueError("fcoefs parameter passed to ERBFilterBank is the wrong size.")

    if x.shape[1] < x.shape[0]:
        x = x.T

    A0 = fcoefs[:, 0]
    A11 = fcoefs[:, 1]
    A12 = fcoefs[:, 2]
    A13 = fcoefs[:, 3]
    A14 = fcoefs[:, 4]
    A2 = fcoefs[:, 5]
    B0 = fcoefs[:, 6]
    B1 = fcoefs[:, 7]
    B2 = fcoefs[:, 8]
    gain = fcoefs[:, 9]

    output = np.zeros((gain.shape[0], x.shape[1]))
    for chan in range(gain.shape[0]):
        y1 = lfilter([A0[chan] / gain[chan], A11[chan] / gain[chan], A2[chan] / gain[chan]],
                     [B0[chan], B1[chan], B2[chan]], x)
        y2 = lfilter([A0[chan], A12[chan], A2[chan]],
                     [B0[chan], B1[chan], B2[chan]], y1)
        y3 = lfilter([A0[chan], A13[chan], A2[chan]],
                     [B0[chan], B1[chan], B2[chan]], y2)
        y4 = lfilter([A0[chan], A14[chan], A2[chan]],
                     [B0[chan], B1[chan], B2[chan]], y3)
        output[chan, :] = y4

    return output


def ERBSpace(lowFreq=100, highFreq=44100/4, N=100):
    # Change the following three parameters if you wish to use a different
    # ERB scale. Must change in MakeERBCoeffs too.
    EarQ = 9.26449  # Glasberg and Moore Parameters
    minBW = 24.7
    order = 1

    # All of the following expressions are derived in Apple TR #35, "An
    # Efficient Implementation of the Patterson-Holdsworth Cochlear
    # Filter Bank." See pages 33-34.
    cfArray = -(EarQ * minBW) + np.exp(
        np.arange(1, N + 1) * (-np.log(highFreq + EarQ * minBW) + np.log(lowFreq + EarQ * minBW)) / N
    ) * (highFreq + EarQ * minBW)

    return cfArray

def MakeERBFilters(fs, numChannels, lowFreq=100):
    T = 1 / fs
    #if len(numChannels) == 1:
    #    cf = ERBSpace(lowFreq, fs / 2, numChannels)
    #else:
    #    cf = numChannels
    #    if cf.shape[1] > cf.shape[0]:
    #        cf = cf.T
    cf = ERBSpace(lowFreq, fs / 2, numChannels)

    # Change the following three parameters if you wish to use a different
    # ERB scale. Must change in ERBSpace too.
    EarQ = 9.26449  # Glasberg and Moore Parameters
    minBW = 24.7
    order = 1

    ERB = ((cf / EarQ) ** order + minBW ** order) ** (1 / order)
    B = 1.019 * 2 * np.pi * ERB

    A0 = T
    A2 = 0
    B0 = 1
    B1 = -2 * np.cos(2 * cf * np.pi * T) / np.exp(B * T)
    B2 = np.exp(-2 * B * T)

    A11 = -(2 * T * np.cos(2 * cf * np.pi * T) / np.exp(B * T) + 2 * np.sqrt(3 + 2 ** 1.5) * T * np.sin(
        2 * cf * np.pi * T) / np.exp(B * T)) / 2
    A12 = -(2 * T * np.cos(2 * cf * np.pi * T) / np.exp(B * T) - 2 * np.sqrt(3 + 2 ** 1.5) * T * np.sin(
        2 * cf * np.pi * T) / np.exp(B * T)) / 2
    A13 = -(2 * T * np.cos(2 * cf * np.pi * T) / np.exp(B * T) + 2 * np.sqrt(3 - 2 ** 1.5) * T * np.sin(
        2 * cf * np.pi * T) / np.exp(B * T)) / 2
    A14 = -(2 * T * np.cos(2 * cf * np.pi * T) / np.exp(B * T) - 2 * np.sqrt(3 - 2 ** 1.5) * T * np.sin(
        2 * cf * np.pi * T) / np.exp(B * T)) / 2

    gain = abs(
        (
                (-2 * np.exp(4j * cf * np.pi * T) * T +
                 2 * np.exp(-(B * T) + 2j * cf * np.pi * T) * T *
                 (np.cos(2 * cf * np.pi * T) - np.sqrt(3 - 2 ** (3 / 2)) *
                  np.sin(2 * cf * np.pi * T)))
                *
                (-2 * np.exp(4j * cf * np.pi * T) * T +
                 2 * np.exp(-(B * T) + 2j * cf * np.pi * T) * T *
                 (np.cos(2 * cf * np.pi * T) + np.sqrt(3 - 2 ** (3 / 2)) *
                  np.sin(2 * cf * np.pi * T)))
                *
                (-2 * np.exp(4j * cf * np.pi * T) * T +
                 2 * np.exp(-(B * T) + 2j * cf * np.pi * T) * T *
                 (np.cos(2 * cf * np.pi * T) - np.sqrt(3 + 2 ** (3 / 2)) *
                  np.sin(2 * cf * np.pi * T)))
                *
                (-2 * np.exp(4j * cf * np.pi * T) * T +
                 2 * np.exp(-(B * T) + 2j * cf * np.pi * T) * T *
                 (np.cos(2 * cf * np.pi * T) + np.sqrt(3 + 2 ** (3 / 2)) *
                  np.sin(2 * cf * np.pi * T)))
        )
        /
        (
                (-2 / np.exp(2 * B * T) - 2 * np.exp(4j * cf * np.pi * T) +
                 2 * (1 + np.exp(4j * cf * np.pi * T)) / np.exp(B * T)) ** 4
        )
    )

    allfilts = np.ones((len(cf), 1))
    A11 = np.expand_dims(A11,axis=1)
    A12 = np.expand_dims(A12,axis=1)
    A13 = np.expand_dims(A13,axis=1)
    A14 = np.expand_dims(A14,axis=1)
    
    B1 = np.expand_dims(B1,axis=1)
    B2 = np.expand_dims(B2,axis=1)
    gain = np.expand_dims(gain,axis=1)
    
    fcoefs = np.concatenate(
        [A0 * allfilts, A11, A12, A13, A14, A2 * allfilts, B0 * allfilts, B1, B2, gain],
        axis=1
    )
    return fcoefs

    if 0:  # Test Code
        A0 = fcoefs[:, 0]
        A11 = fcoefs[:, 1]
        A12 = fcoefs[:, 2]
        A13 = fcoefs[:, 3]
        A14 = fcoefs[:, 4]
        A2 = fcoefs[:, 5]
        B0 = fcoefs[:, 6]
        B1 = fcoefs[:, 7]
        B2 = fcoefs[:, 8]
        gain = fcoefs[:, 9]
        chan = 1
        x = np.concatenate(([1], np.zeros(511)))
        y1 = np.convolve([A0[chan] / gain[chan], A11[chan] / gain[chan], A2[chan] / gain[chan]], [B0[chan], B1[chan], B2[chan]], x)
        y2 = np.convolve([A0[chan], A12[chan], A2[chan]], [B0[chan], B1[chan], B2[chan]], y1)
        y3 = np.convolve([A0[chan], A13[chan], A2[chan]], [B0[chan], B1[chan], B2[chan]], y2)
        y4 = np.convolve([A0[chan], A14[chan], A2[chan]], [B0[chan], B1[chan], B2[chan]], y3)
        freqScale = np.arange(len(x)) * (fs / len(x))
        resp = 20 * np.log10(np.abs(np.fft.fft(y4)))
        semilogx(freqScale[:256], resp[:256])
        
def AKerbILD_short_p1(P, f_lim, fs):
    x_l = np.squeeze(P[:,0,:])
    x_l = np.transpose(x_l,(1,0))
    x_r = np.squeeze(P[:,1,:])
    x_r = np.transpose(x_r,(1,0))
    
    # get length and channels
    N = x_l.shape[0]  # ir length
    c = x_l.shape[1]  # ir num. of channels
    

    # ---------------------------------------- 1. get FIR filter in freq domain
    # get ERB filter coefficients
    n      = round(21.4*np.log10(0.004367*f_lim[1]+1)-21.4*np.log10(0.004367*f_lim[0]+1))    # number of filters
    f_c    = ERBSpace(f_lim[0], fs/2, n)                                               # center frequencies
    fcoefs = MakeERBFilters(fs,n,f_lim[0])  
    
    # Sort with ascending frequency
    fcoefs = np.flipud(fcoefs)
    f_c = np.flipud(f_c)

    # Get rid of filters above upper f_lim
    n = np.where(f_c <= f_lim[1])[0][-1] + 1
    f_c = f_c[:n]
    fcoefs = fcoefs[:n, :]

    # Estimate needed filter length (4 times the largest cycle)
    Nmax = int(2 ** np.ceil(np.log2(1 / f_c[0] * 4 * fs)))

    # Get filter impulse responses
    C = np.zeros((Nmax,1))
    C[0] = 1
    C = np.transpose(ERBFilterBank(np.transpose(C),fcoefs));

    # Make minimum phase
    # C = AKphaseManipulation(C, fs, 'min', 4, 0)  # Uncomment if necessary

    # Filter x with filter bank and calculate error in each band
    ILD = np.zeros((n, c))
    C_db_tmp = 10 * np.log10(np.abs(C))
    C_db = np.zeros(C_db_tmp.shape)
    
    return Nmax,f_c,n,C
    


