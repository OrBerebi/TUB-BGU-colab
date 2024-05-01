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
        p_hat_t = torch.fft.irfft(p, n=filt_samp, dim=1)  # space x time x left/right
        p_hat_t = torch.roll(p_hat_t, shifts=filt_samp // 2, dims=1)
        return p_hat_t
        
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
        e_nmse_r =  torch.linalg.vector_norm(p_ref_r - p_r, dim=0)
        
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

def clc_e_mag_v2(p_ref,p_hat,M_dB,norm_flag=True):

    #calc the absolute mean diffrence
    e_mag_space_freq_lr = torch.abs(torch.abs(p_ref) - torch.abs(p_hat))

    # mask the only posirive gradient bins
    if M_dB.numel() > 0: 
        #The tensor is not empty
        M_linear = 10 ** (M_dB / 10)
        e_mag_space_freq_lr = e_mag_space_freq_lr * M_linear 
    

    # norm and avarage over all directions ears and freqncies
    if norm_flag:
        e_mag_freq_lr = torch.norm(e_mag_space_freq_lr, dim=0).squeeze() / torch.norm(torch.abs(p_ref), dim=0).squeeze() # || || over directions
    else:
        e_mag_freq_lr = torch.norm(e_mag_space_freq_lr, dim=0).squeeze() # || || over directions
        
    
    e_mag_freq_lr_db = 20*torch.log10(e_mag_freq_lr)                 # calc dB error
    e_mag_freq = torch.mean(e_mag_freq_lr_db,dim = 1).squeeze()      # E[] over ears
    
    return e_mag_freq


