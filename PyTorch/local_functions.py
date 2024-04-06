import torch
import numpy as np
from scipy.fft import rfft,irfft
from scipy.signal import lfilter

def spatial_interpolation(Hnm,Y,is_time=False):
    filt_samp  = int((Hnm.shape[1] -1)* 2)
    p_left     = torch.matmul(Y, Hnm[:,:,0])
    p_right    = torch.matmul(Y, Hnm[:,:,1])
    p          = torch.stack((p_left,p_right),2)

    if is_time:
        p_t = np.fft.irfft(p, n = filt_samp,axis = 1) # space x time x left/right
        p_t = np.roll(p_t, shift=filt_samp // 2, axis = 1)
        return(p_t)
        
    return(p)

def clc_e_nmse(p_ref,p_hat,norm_flag=True):
    p_l = p_hat[:,:,0]
    p_r = p_hat[:,:,1]
    
    p_ref_l = p_ref[:, :, 0]
    p_ref_r = p_ref[:, :, 1]
    
    if norm_flag:
        e_nmse_l =  torch.linalg.vector_norm(p_ref_l - p_l, dim=0) / torch.linalg.vector_norm(p_ref_l, dim=0)
        e_nmse_r =  torch.linalg.vector_norm(p_ref_r - p_r, dim=0) / torch.linalg.vector_norm(p_ref_r, dim=0)
    else:
        e_nmse_l =  torch.linalg.vector_norm(p_ref_l - p_l, dim=0)
        e_nmse_r = torch.linalg.vector_norm(p_ref_r - p_r, dim=0)
        
    e_nmse_norm = 20 * torch.log10((e_nmse_l + e_nmse_r) / 2)

    return e_nmse_norm

def clc_e_mag(p_ref,p_hat,norm_flag=True):
    p_l = p_hat[:,:,0]
    p_r = p_hat[:,:,1]
    
    p_ref_l = p_ref[:, :, 0]
    p_ref_r = p_ref[:, :, 1]

    if norm_flag:
        e_mag_l = torch.norm(torch.abs(p_ref_l) - torch.abs(p_l), dim=0) / torch.norm(torch.abs(p_ref_l), dim=0)
        e_mag_r = torch.norm(torch.abs(p_ref_r) - torch.abs(p_r), dim=0) / torch.norm(torch.abs(p_ref_r), dim=0)
    else:
        e_mag_l = torch.norm(torch.abs(p_ref_l) - torch.abs(p_l), dim=0)
        e_mag_r = torch.norm(torch.abs(p_ref_r) - torch.abs(p_r), dim=0) 

    e_mag_norm = 20 * torch.log10((e_mag_l + e_mag_r) / 2)


    return e_mag_norm


