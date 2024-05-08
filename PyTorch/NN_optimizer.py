
import os
import sys
import time
import platform
import torch
import torchaudio

sys.path.append('/Users/orberebi/Documents/GitHub/TUB-BGU-colab/PyTorch/')
import scipy.io
import numpy as np
np.set_printoptions(precision=3)
np.set_printoptions(threshold=4)
import torch
import local_functions as f
import AK_ILD_pytorch as ak
import ITD_estimator as itd_est


from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import pyfar as pf
import sofar as sf

from torchinfo import summary

import baumgartner2014 as bam
import pickle




class NN_v1(torch.nn.Module):
    def __init__(self,M,freq_bins):
        super().__init__()
        self.linear1        = torch.nn.Linear(M*2,M*2)
        self.linear1.weight = torch.nn.Parameter(torch.eye(M*2,dtype=torch.complex128), requires_grad=True)
        self.linear1.bias   = torch.nn.Parameter(torch.zeros((freq_bins,M*2),dtype=torch.complex128), requires_grad=True)
        
        self.linear2        = torch.nn.Linear(freq_bins,freq_bins)
        self.linear2.weight = torch.nn.Parameter(torch.eye(freq_bins,dtype=torch.complex128), requires_grad=True)
        self.linear2.bias   = torch.nn.Parameter(torch.zeros((M*2,freq_bins),dtype=torch.complex128), requires_grad=True)
        
        self.linear3        = torch.nn.Linear(M*2,M*2)
        self.linear3.weight = torch.nn.Parameter(torch.eye(M*2,dtype=torch.complex128), requires_grad=True)
        self.linear3.bias   = torch.nn.Parameter(torch.zeros((freq_bins,M*2),dtype=torch.complex128), requires_grad=True)
        
        self.linear4        = torch.nn.Linear(freq_bins,freq_bins)
        self.linear4.weight = torch.nn.Parameter(torch.eye(freq_bins,dtype=torch.complex128), requires_grad=True)
        self.linear4.bias   = torch.nn.Parameter(torch.zeros((M*2,freq_bins),dtype=torch.complex128), requires_grad=True)
                
        self.act_func = complex_tanh_layer()
        
        self.M = M
        self.freq_bins = freq_bins
        

    def forward(self, x):
        
        x  = torch.reshape(x,(self.freq_bins, 2*self.M))
        
        h1   = self.linear1(x)
        h1   = self.act_func(h1)
        
        h1_t = torch.transpose(h1, 0, 1)
        
        h2   = self.linear2(h1_t)
        h2   = self.act_func(h2)
        h2_t = torch.transpose(h2, 0, 1)
  
        h1   = self.linear3(h2_t)
        h1_t = torch.transpose(h1, 0, 1)
        h2   = self.linear4(h1_t)
        h2_t = torch.transpose(h2, 0, 1)
        
        #y = torch.reshape(h2_t,(self.freq_bins, self.M, 2))
        y = torch.reshape(h2_t,(self.M, self.freq_bins, 2))
        
        return y



class NN_v2_simp(torch.nn.Module):
    def __init__(self,M,freq_bins):
        super().__init__()
        self.linear1        = torch.nn.Linear(M*2,M*2)
        self.linear1.weight = torch.nn.Parameter(torch.eye(M*2,dtype=torch.complex128), requires_grad=True)
        self.linear1.bias   = torch.nn.Parameter(torch.zeros((freq_bins,M*2),dtype=torch.complex128), requires_grad=True)
        
        self.linear2        = torch.nn.Linear(freq_bins,freq_bins)
        self.linear2.weight = torch.nn.Parameter(torch.eye(freq_bins,dtype=torch.complex128), requires_grad=True)
        self.linear2.bias   = torch.nn.Parameter(torch.zeros((M*2,freq_bins),dtype=torch.complex128), requires_grad=True)
        
        self.M = M
        self.freq_bins = freq_bins
        

    def forward(self, x):
        
        x  = torch.reshape(x,(self.freq_bins, 2*self.M))
        
        h1   = self.linear1(x)        
        h1_t = torch.transpose(h1, 0, 1)
        
        h2   = self.linear2(h1_t)
        h2_t = torch.transpose(h2, 0, 1)
        
        y = torch.reshape(h2_t,(self.M, self.freq_bins, 2))
        
        return y

class complex_tanh_layer(torch.nn.Module):
    def __init__(self):
        """
        In the constructor we instantiate four parameters and assign them as
        member parameters.
        """
        super().__init__()

    def forward(self, input):
        """
        In the forward function we accept a Tensor of input data and we must return
        a Tensor of output data. We can use Modules defined in the constructor as
        well as arbitrary operators on Tensors.
        """
        return torch.tanh(input.real).type(torch.complex128)+1j*torch.tanh(input.imag).type(torch.complex128)

def calc_e_enrgy(p_hat,p_ref):
    enrgy_ref = torch.squeeze(torch.norm(torch.abs(p_ref[:,:,0]), dim=1)) / torch.squeeze(torch.norm(torch.abs(p_ref[:,:,1]), dim=1))
    enrgy_hat = torch.squeeze(torch.norm(torch.abs(p_hat[:,:,0]), dim=1)) / torch.squeeze(torch.norm(torch.abs(p_hat[:,:,1]), dim=1))
    
    e_enrgy = torch.squeeze(torch.norm(torch.abs(enrgy_ref-enrgy_hat), dim=0))
    return e_enrgy

def clc_e_mag_diff(p_hat):
    p_l = torch.diff(torch.abs(p_hat[:,:,0]),n=1, dim=1)
    p_r = torch.diff(torch.abs(p_hat[:,:,1]),n=1, dim=1)
    e_mag_diff_l = torch.norm(torch.abs(p_l), dim=0)
    e_mag_diff_r = torch.norm(torch.abs(p_r), dim=0) 

    e_mag_norm = ((e_mag_diff_l + e_mag_diff_r) / 2)


    return e_mag_norm


def loss_function_v1(Hnm,Y_lebedev,Y_az,p_ref,ILD_ref, AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft,lambda_vec,norm_flag):
    # Calc the low-order HRTF in the Omega directions
    p_lebedev = f.spatial_interpolation(Hnm,Y_lebedev,False)
    
    # Calculate low order p in time (to evaluate it's ILD)
    p_az_t  = f.spatial_interpolation(Hnm,Y_az,True)
    
    # Calc the NMSE and Magnitude error for each set of Hnm
    e_nmse = f.clc_e_nmse(p_ref,p_lebedev,norm_flag)
    e_mag  = f.clc_e_mag(p_ref,p_lebedev,norm_flag)
    
    # Calculate the ILD error for each low-order representation
    ILD = ak.clc_ILD(p_az_t, AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)
    e_ILD    = torch.mean(torch.abs(ILD_ref - ILD),dim=0) #avarage over frequencies bands

    # Get error band indexs
    f_vec   = np.arange(0,e_mag.shape[0])*((fs/2)/(e_mag.shape[0]-1))
    idx_min = (np.abs(f_vec - f_band[0])).argmin()
    idx_max = (np.abs(f_vec - f_band[1])).argmin()
    
    idx_max_mag = (np.abs(f_vec - fs/2)).argmin()
    idx_min_mag = (np.abs(f_vec - 2e3)).argmin()
    
    # Calc error over frequencies and over directions
    lambda_0 = lambda_vec[0]
    lambda_1 = lambda_vec[1]
    lambda_2 = lambda_vec[2]

    e_total  = lambda_0*torch.norm(e_ILD,dim=0) + lambda_1*torch.norm(e_nmse,dim=0) + lambda_2*torch.norm(e_mag[idx_min_mag:idx_max_mag],dim=0)

    return e_total, e_ILD, e_nmse, e_mag, ILD


def loss_function_v2(Hnm,Y_lebedev,Y_az,p_ref,p_ref_az_t ,ILD_ref, AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft,lambda_vec,norm_flag):
    # Calc the low-order HRTF in the Omega directions
    p_lebedev = f.spatial_interpolation(Hnm,Y_lebedev,False)
    
    # Calculate low order p in time (to evaluate it's ILD)
    p_az_t  = f.spatial_interpolation(Hnm,Y_az,True)
    
    # Calc the NMSE and Magnitude error for each set of Hnm
    e_nmse = f.clc_e_nmse(p_ref,p_lebedev,norm_flag)
    e_mag  = f.clc_e_mag(p_ref,p_lebedev,norm_flag)
    
    # Calculate the ILD error for each low-order representation
    ILD      = ak.clc_ILD(p_az_t, AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)
    e_ILD    = torch.mean(torch.abs(ILD_ref - ILD),dim=0) #avarage over frequencies bands

    # Calculate the Enargy error for each low-order representation
    e_mag_diff = clc_e_mag_diff(p_lebedev)
    
    # Get error band indexs
    f_vec   = np.arange(0,e_mag.shape[0])*((fs/2)/(e_mag.shape[0]-1))
    idx_min = (np.abs(f_vec - f_band[0])).argmin()
    idx_max = (np.abs(f_vec - f_band[1])).argmin()
    
    idx_max_mag = (np.abs(f_vec - fs/2)).argmin()
    idx_min_mag = (np.abs(f_vec - 2e3)).argmin()
    
    # Calc error over frequencies and over directions
    lambda_0 = lambda_vec[0]
    lambda_1 = lambda_vec[1]
    lambda_2 = lambda_vec[2]
    lambda_3 = lambda_vec[3]

    e_total  =  lambda_0*torch.norm(e_ILD,dim=0) + lambda_1*torch.norm(e_nmse,dim=0) + lambda_2*torch.norm(e_mag[idx_min_mag:idx_max_mag],dim=0) + lambda_3*torch.norm(e_mag_diff[idx_min_mag:idx_max_mag],dim=0) 

    return e_total, e_ILD, e_nmse, e_mag, ILD, torch.norm(e_mag_diff[idx_min_mag:idx_max_mag],dim=0)








def loss_function_v3(Hnm,Y_lebedev,Y_az,p_ref,p_ref_az_t ,ILD_ref, AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft,lambda_vec,norm_flag):
    # Calc the low-order HRTF in the Omega directions
    p_lebedev = f.spatial_interpolation(Hnm,Y_lebedev,False)
    
    # Calculate low order p in time (to evaluate it's ILD)
    p_az_t  = f.spatial_interpolation(Hnm,Y_az,True)
    
    # Calc the NMSE and Magnitude error for each set of Hnm
    e_nmse = f.clc_e_nmse(p_ref,p_lebedev,norm_flag)
    e_mag  = f.clc_e_mag(p_ref,p_lebedev,norm_flag)
    
    # Calculate the ILD error for each low-order representation
    ILD      = ak.clc_ILD(p_az_t, AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)
    e_ILD    = torch.mean(torch.abs(ILD_ref - ILD),dim=0) #avarage over frequencies bands

    # Calculate colorization error
    e_color_fc_lr_space = ak.calc_color(p_ref_az_t,p_az_t,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)
    e_color_fc_space    = torch.mean(e_color_fc_lr_space,dim=2)
    e_color_fc          = torch.mean(e_color_fc_space,dim=1)
    e_color             = torch.abs(1 - torch.abs(torch.mean(e_color_fc,dim=0))) # good if close to zero
    
    

    # Calculate the Enargy error for each low-order representation
    e_mag_diff = clc_e_mag_diff(p_lebedev)
    
    # Get error band indexs
    f_vec   = np.arange(0,e_mag.shape[0])*((fs/2)/(e_mag.shape[0]-1))
    idx_min = (np.abs(f_vec - f_band[0])).argmin()
    idx_max = (np.abs(f_vec - f_band[1])).argmin()
    
    idx_max_mag = (np.abs(f_vec - fs/2)).argmin()
    idx_min_mag = (np.abs(f_vec - 2e3)).argmin()
    
    # Calc error over frequencies and over directions
    lambda_0 = lambda_vec[0]
    lambda_1 = lambda_vec[1]
    lambda_2 = lambda_vec[2]
    lambda_3 = lambda_vec[3]
    lambda_4 = lambda_vec[4]

    e_total  =  lambda_0*torch.norm(e_ILD,dim=0) + lambda_1*torch.norm(e_nmse,dim=0) + lambda_2*torch.norm(e_mag[idx_min_mag:idx_max_mag],dim=0) + lambda_3*torch.norm(e_mag_diff[idx_min_mag:idx_max_mag],dim=0) + lambda_4*e_color
    
    e_mag_diff = torch.norm(e_mag_diff[idx_min_mag:idx_max_mag],dim=0)

    output_dict = {
        "e_total": e_total,
        "e_ILD": e_ILD,
        "e_nmse": e_nmse,
        "e_mag": e_mag,
        "ILD": ILD,
        "e_mag_diff": e_mag_diff,
        "e_color": e_color
    }

    return output_dict

def loss_function_v4(Hnm,Y_lebedev,Y_az,p_ref,p_ref_az_t ,ILD_ref,ITD_ref,GD_ref,pos_grad_ref_fc,pos_grad_ref_f, AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft,lambda_vec,save_path_sofa,SourcePosition,flag_s2,norm_flag):
    # Calc the low-order HRTF in the Omega directions
    p_lebedev = f.spatial_interpolation(Hnm,Y_lebedev,False)
    
    # Calculate low order p in time (to evaluate it's ILD)
    p_az_t  = f.spatial_interpolation(Hnm,Y_az,True)

    if flag_s2:
        p_lebedev_t  = f.spatial_interpolation(Hnm,Y_lebedev,True)

    #==============================================================================
    # Calc the NMSE error for each set of Hnm
    #==============================================================================
    if lambda_vec[1] == 0:
        e_nmse = torch.zeros(Hnm.shape[1], requires_grad=True)
    else:
        e_nmse = f.clc_e_nmse(p_ref,p_lebedev,norm_flag)

    #==============================================================================
    # Calc the Magnitude error for each set of Hnm
    #==============================================================================
    if lambda_vec[7] == 0:
        e_mag_freq_weighted = torch.zeros(Hnm.shape[1], requires_grad=True)
    else:
        e_mag_freq_weighted = f.clc_e_mag_v2(p_ref,p_lebedev,pos_grad_ref_f,norm_flag) # [space x freq x ears]
    
    if lambda_vec[2] == 0:
        e_mag_freq = torch.zeros(Hnm.shape[1], requires_grad=True)
    else:
        e_mag_freq  = f.clc_e_mag_v2(p_ref,p_lebedev,torch.tensor([]),norm_flag) # [space x freq x ears]


    #==============================================================================
    # Calculate the ILD error for each low-order representation
    #==============================================================================
    if flag_s2:
        if lambda_vec[0] == 0:
            ILD = torch.zeros(AK_f_c.shape[0], p_lebedev_t.shape[0] , requires_grad=True)
        else:
            ILD      = ak.clc_ILD(p_lebedev_t, AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)
    else:
        if lambda_vec[0] == 0:
            ILD = torch.zeros(AK_f_c.shape[0], p_az_t.shape[0] , requires_grad=True)
        else:
            ILD      = ak.clc_ILD(p_az_t, AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)
    
    e_ILD    = torch.mean(torch.abs(ILD_ref - ILD),dim=0) #avarage over frequencies bands

    #==============================================================================
    # Calculate the ITD error for low-order representation
    #==============================================================================
    if lambda_vec[5] == 0:
        ITD = torch.zeros(p_az_t.shape[0] , requires_grad=True)
        GD  = torch.zeros(GD_ref.shape , requires_grad=True)
    else:
        ITD   = itd_est.get_ITD_group_delay(p_az_t.float(),fs, 3e3,False) #returns the micro sec ITD usinng the onset method
        GD   = itd_est.get_muilti_ch_group_delay(p_az_t.float(),fs, 3e3,False) # returns the group delay in [freq x space x ears]
        
    e_GD = torch.mean(torch.abs(GD_ref - GD),dim= 1 ).squeeze() # avarage over directions
    e_GD = torch.mean(e_GD,dim= 1).squeeze() # avarage over ears
    e_ITD = torch.norm(e_GD,dim= 0) # norm over freqncies

    #e_ITD = torch.mean(torch.abs(ITD_ref - ITD)) # avarage over freqncies
    
    #print("group delay error value: ", e_GD)
    #print("ITD error value: ", e_ITD)
    
    #==============================================================================
    # Calculate colorization error
    #==============================================================================
    if lambda_vec[4] == 0:
        e_color_fc_lr_space = torch.zeros(AK_f_c.shape[0], p_lebedev.shape[0] ,p_lebedev.shape[2]  , requires_grad=True)
    else:
        e_color_fc_lr_space = ak.calc_color(p_ref,p_lebedev,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)

    if lambda_vec[6] == 0:
        e_color_fc_lr_space_pos_grad_weighted = torch.zeros(AK_f_c.shape[0], p_lebedev.shape[0] ,p_lebedev.shape[2]  , requires_grad=True)
    else:
        e_color_fc_lr_space_pos_grad_weighted = ak.calc_color(p_ref,p_lebedev,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)
        pos_grad_ref_fc_linear = 10 ** (pos_grad_ref_fc / 20) # the mask is given in dB turn to linear for multiplication
        e_color_fc_lr_space_pos_grad_weighted = e_color_fc_lr_space_pos_grad_weighted * pos_grad_ref_fc_linear # mask the "important" directions and frequencies
        e_color_fc_lr_space_pos_grad_weighted = (10*torch.log10(torch.abs(e_color_fc_lr_space_pos_grad_weighted)))**2 # zero is good


    
    e_color_fc_lr_space_out = e_color_fc_lr_space

    e_color_fc_lr_space = (10*torch.log10(e_color_fc_lr_space))**2 # zero is good
    
    e_color_fc_lr_weighted  = torch.mean(e_color_fc_lr_space_pos_grad_weighted,dim=1).squeeze() # E[.] over space
    e_color_fc_weighted     = torch.mean(e_color_fc_lr_weighted,dim=1).squeeze()                # E[.] over ears
    e_color_weighted        = torch.mean(e_color_fc_weighted,dim=0).squeeze()                   # E[.] over freqncies
    
    e_color_fc_lr  = torch.mean(e_color_fc_lr_space,dim=1).squeeze() # E[.] over space
    e_color_fc     = torch.mean(e_color_fc_lr,dim=1).squeeze()       # E[.] over ears
    e_color        = torch.mean(e_color_fc,dim=0).squeeze()          # E[.] over freqncies
    
    
    #==============================================================================
    # Calculate the derivitive magnitude, to smooth results
    #==============================================================================
    if lambda_vec[3] == 0:
        e_mag_diff = torch.zeros( Hnm.shape[1] - 1  , requires_grad=True)
    else:
        e_mag_diff = clc_e_mag_diff(p_lebedev)
    
    
    
    # Get error band indexs
    f_vec   = np.arange(0,e_nmse.shape[0])*((fs/2)/(e_nmse.shape[0]-1))
    idx_min = (np.abs(f_vec - f_band[0])).argmin()
    idx_max = (np.abs(f_vec - f_band[1])).argmin()
    
    idx_max_mag = (np.abs(f_vec - fs/2)).argmin()
    idx_min_mag = (np.abs(f_vec - 2e3)).argmin()
    
    # Calc error over frequencies and over directions
    lambda_0 = lambda_vec[0]
    lambda_1 = lambda_vec[1]
    lambda_2 = lambda_vec[2]
    lambda_3 = lambda_vec[3]
    lambda_4 = lambda_vec[4]
    lambda_5 = lambda_vec[5]
    lambda_6 = lambda_vec[6]
    lambda_7 = lambda_vec[7]

    # Save results as sofa for the baumgartner2014 model
    HRIR_name = "/iMagLS.sofa"
    HRIR_ref_name = "/ref.sofa"
    Hnm_out = Hnm.detach().cpu()
    save_as_sofa(Hnm_out,Y_lebedev,SourcePosition,int(fs),save_path_sofa,HRIR_name)
    do_dtf = True
    shutup = True
    sofa_path_target = save_path_sofa + HRIR_name
    sofa_path_template = save_path_sofa + HRIR_ref_name
    sig_path = []
    out = bam.baumgartner2014(sofa_path_target,sofa_path_template,sig_path,shutup,do_dtf)
    local_pe, local_qe, circ_mean, circ_var, circ_std = bam.circular_stats_from_pdf(out['pdf'], out['rang'], out['tang'],shutup)
    local_qe = np.mean(local_qe)
    local_pe = np.mean(local_pe)
    

    # Define the individual terms
    term1 = lambda_0 * torch.norm(e_ILD, dim=0)
    term2 = lambda_1 * torch.norm(e_nmse, dim=0)
    term3 = lambda_2 * torch.norm(e_mag_freq[idx_min_mag:idx_max_mag], dim=0)
    term4 = lambda_3 * torch.norm(e_mag_diff[idx_min_mag:idx_max_mag], dim=0)
    term5 = lambda_4 * e_color
    term6 = lambda_5 * e_ITD
    term7 = lambda_6 * e_color_weighted
    term8 = lambda_7 * torch.norm(e_mag_freq_weighted[idx_min_mag:idx_max_mag], dim=0)
    
    # Sum up the terms
    e_total = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8
    
    e_mag_diff = torch.norm(e_mag_diff[idx_min_mag:idx_max_mag],dim=0)

    output_dict = {
        "e_total": e_total,
        "e_ILD": e_ILD,
        "e_nmse": e_nmse,
        "e_mag": e_mag_freq,
        "ILD": ILD,
        "e_mag_diff": e_mag_diff,
        "e_color": e_color,
        "e_color_fc_lr_space": e_color_fc_lr_space_out,
        "e_ITD": e_ITD,
        "qe": local_qe,
        "pe": local_pe,
        "e_color_weighted": e_color_weighted,
        "e_mag_weighted": e_mag_freq_weighted
    }


    return output_dict


def calc_positive_grad(Hnm,Y,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft,shutup):
    p_space_f_lr = f.spatial_interpolation(Hnm,Y,True)
    x = ak.calc_gammatone_spectrum(p_space_f_lr,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)

    x = 10 * torch.log10(x)
    x_pos_grad = torch.zeros(x.shape[0], x.shape[1], x.shape[2], dtype=x.dtype, device=x.device)
    
    
    # Calculate the gradient
    for f_idx in range(x.shape[0] - 1):
        x_pos_grad[f_idx, :, :] = x[f_idx + 1, :, :] - x[f_idx, :, :]
    
    # Take only the positive values
    x_pos_grad[x_pos_grad < 0] = 0

    if not(shutup):
        # Display using imshow
        plt.imshow(x_pos_grad[:,::5, 0].detach().cpu().numpy(), cmap='bone', origin='lower')  # Assuming the first dimension represents frequency
        plt.colorbar()  # Add colorbar
        plt.title("Positive Gradient of ERB spectrum")
        plt.xlabel("angle indices")
        plt.ylabel("ERB filter idex")
        plt.show()

    
    return x_pos_grad

def do_rms_pf(sig):
    # Averaging over time (RMS)
    # Calculate root mean square (rms) along the last axis
    rms_values = np.sqrt(np.mean(np.square(sig.time), axis=-1))
    # Squeeze the array to remove axes with length 1
    rms_values_squeezed = np.squeeze(rms_values)
    # Convert to dB scale
    sig_rms_dB = 20 * np.log10(rms_values_squeezed)
    return sig_rms_dB
    
def calc_positive_grad_pyfar(Hnm,Y,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft,shutup):
    p_space_f_lr = f.spatial_interpolation(Hnm,Y,True)

    x = p_space_f_lr.permute(2,0,1)
    x = pf.Signal(x,fs)
    Gammatones     = pf.dsp.filter.GammatoneBands(freq_range=[20, 20e3],sampling_rate=fs) # Create the filter bank
    x         = Gammatones.process(x)[0] # Apply the filter bank (T)
    x   = do_rms_pf(x)

    AK_f_c = Gammatones.frequencies
    
    #x = ak.calc_gammatone_spectrum(p_space_f_lr,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)
    #x = 10 * torch.log10(x)
    x_pos_grad = np.zeros([x.shape[0],x.shape[1],x.shape[2]])
    # Calculate the gradient
    for f_idx in range(x.shape[0]-1):
        x_pos_grad[f_idx,:,:]   = x[f_idx+1,:,:] - x[f_idx,:,:]
    # Take only the positive values
    x_pos_grad[x_pos_grad < 0] = 0



    x_pos_grad =  np.transpose(x_pos_grad, (0, 2, 1))

    poo = np.transpose(x_pos_grad, (2, 1, 0))
    x_pos_grad_pf = pf.FrequencyData(poo, AK_f_c)

    x_pos_grad_pf_left  = x_pos_grad_pf[0,:,:]
    x_pos_grad_pf_right = x_pos_grad_pf[1,:,:]

    #out_freq_size = int(p_space_f_lr.shape[1]/2 + 1)
    out_freq_size = int(p_space_f_lr.shape[1])
    
    interpolator = pf.dsp.InterpolateSpectrum(x_pos_grad_pf_left, method = 'magnitude', kind = ('nearest', 'linear', 'nearest'),fscale= 'log')
    x_pos_grad_interpulated = interpolator(out_freq_size, fs,show=False)
    left_ear_inteprulated = x_pos_grad_interpulated.freq.real
    
    interpolator = pf.dsp.InterpolateSpectrum(x_pos_grad_pf_right, method = 'magnitude', kind = ('nearest', 'linear', 'nearest'),fscale= 'log')
    x_pos_grad_interpulated = interpolator(out_freq_size, fs,show=False)
    right_ear_inteprulated = x_pos_grad_interpulated.freq.real

    # Reshape the arrays to add a new axis
    #left_ear_inteprulated = left_ear_inteprulated[:, :, np.newaxis]
    #right_ear_inteprulated = right_ear_inteprulated[:, :, np.newaxis]
    
    # Stack the arrays along the third dimension
    inteprulated = np.stack((left_ear_inteprulated, right_ear_inteprulated), axis=2)
    inteprulated = np.transpose(inteprulated, (1, 0, 2))


    f_vec = np.linspace(0,fs/2,inteprulated.shape[0])
    
    idx_max_mag = (np.abs(f_vec - AK_f_c[-1])).argmin()
    idx_min_mag = (np.abs(f_vec - AK_f_c[0])).argmin()


    #inteprulated[0:idx_min_mag,:,:] = 0
    #inteprulated[idx_max_mag:-1,:,:] = 0

    #x_pos_grad = inteprulated




    if not(shutup):
        # Display using imshow
        plt.imshow(inteprulated[:,:, 0], cmap='bone', origin='lower')  # Assuming the first dimension represents frequency
        plt.colorbar()  # Add colorbar
        plt.title("Positive Gradient of bin spectrum")
        plt.xlabel("angle indices")
        plt.ylabel("nfft filter index (positive)")
        plt.show()

        plt.figure
        plt.imshow(x_pos_grad[:,::5, 0], cmap='bone', origin='lower')  # Assuming the first dimension represents frequency
        plt.colorbar()  # Add colorbar
        plt.title("Positive Gradient of ERB spectrum")
        plt.xlabel("angle indices")
        plt.ylabel("ERB filter idex")
        plt.show()


   
    

    

    inteprulated_torch = torch.tensor(inteprulated)
    inteprulated_torch = inteprulated_torch.permute(1,0,2) # [space x freq x ears]
    return inteprulated_torch
    


def calc_positive_grad_bins(Hnm,Y,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft,shutup):
    p_space_f_lr = f.spatial_interpolation(Hnm,Y,True)
    #x = ak.calc_gammatone_spectrum(p_space_f_lr,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)
    x = p_space_f_lr

    x = 20 * torch.log10(torch.abs(x))
    x = x.permute(1, 0, 2)
    print(x.shape)
    x_pos_grad = torch.zeros(x.shape[0], x.shape[1], x.shape[2], dtype=x.dtype, device=x.device)
    
    
    # Calculate the gradient
    for f_idx in range(x.shape[0] - 1):
        x_pos_grad[f_idx, :, :] = x[f_idx + 1, :, :] - x[f_idx, :, :]
    
    # Take only the positive values
    x_pos_grad[x_pos_grad < 0] = 0

    if not(shutup):
        # Display using imshow
        plt.imshow(x_pos_grad[:,::1, 0].detach().cpu().numpy(), cmap='bone', origin='lower')  # Assuming the first dimension represents frequency
        plt.colorbar()  # Add colorbar
        plt.title("Positive Gradient on frequncy spectrum")
        plt.xlabel("angle indices")
        plt.ylabel("nfft bins idex")
        plt.show()

    
    return x_pos_grad

def calc_positive_grad_interpulate(Hnm,Y,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft,shutup):
    p_space_f_lr = f.spatial_interpolation(Hnm,Y,True)
    x = ak.calc_gammatone_spectrum(p_space_f_lr,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)

    x = 10 * torch.log10(x)

    print(p_space_f_lr.shape)
    # Calculate the interpolation factor
    #scale_factor = p_space_f_lr.size(1) / x.size(0)

    # mini-batch x channels x [optional depth] x [optional height] x width.
    print("before: ", x.shape)

    out_freq_size = int(p_space_f_lr.shape[1]/2 + 1)
    x = x.permute(2, 1, 0)
    # Perform interpolation

    x = torch.nn.functional.interpolate(x, size=out_freq_size, mode='linear', align_corners=False)
    x = x.permute(2, 1, 0)

    print("after: ", x.shape)
    x_pos_grad = torch.zeros(x.shape[0], x.shape[1], x.shape[2], dtype=x.dtype, device=x.device)
    
    
    # Calculate the gradient
    for f_idx in range(x.shape[0] - 1):
        x_pos_grad[f_idx, :, :] = x[f_idx + 1, :, :] - x[f_idx, :, :]
    
    # Take only the positive values
    x_pos_grad[x_pos_grad < 0] = 0

    f_vec = np.linspace(0,fs/2,out_freq_size)
    plt.figure
    plt.plot(AK_f_c,linestyle='dashed', marker='o')
    plt.plot(f_vec,linestyle='dashed', marker='o')
    plt.show()
    if not(shutup):
        # Display using imshow
        plt.imshow(x_pos_grad[:,::1, 0].detach().cpu().numpy(), cmap='bone', origin='lower')  # Assuming the first dimension represents frequency
        plt.colorbar()  # Add colorbar
        plt.title("Positive Gradient on frequncy spectrum")
        plt.xlabel("angle indices")
        plt.ylabel("nfft bins idex")
        plt.show()

    
    return x_pos_grad


# Assuming your matrix is named 'input_matrix'
input_matrix = torch.randn(22, 361)  # Replace this with your actual data

# Define the desired output size
output_size = (512, 361)



# Now, interpolated_matrix has dimensions [512, 361]




def do_rms(sig):
    # Averaging over time (RMS)
    # Calculate root mean square (rms) along the last axis
    rms_values = np.sqrt(np.mean(np.square(sig.time), axis=-1))
    # Squeeze the array to remove axes with length 1
    rms_values_squeezed = np.squeeze(rms_values)
    # Convert to dB scale
    sig_rms_dB = 20 * np.log10(rms_values_squeezed)
    return sig_rms_dB

def binaural_reproduction(Hnm,Y,sig):
    
    HRIR = f.spatial_interpolation(Hnm,Y,True)
    HRIR = torch.permute(HRIR,[0,2,1])

    sig = torch.unsqueeze(sig,0)
    sig = sig.expand(HRIR.shape[0], sig.shape[1], sig.shape[2])
    sig = sig.repeat(1, HRIR.shape[1], 1)

    hrir_example_l = torch.squeeze(HRIR[120,0,:])
    hrir_example_r = torch.squeeze(HRIR[120,1,:])

    
    plt.plot(hrir_example_l)
    plt.plot(hrir_example_r)
    plt.ylim((-1, 1))
    plt.show()
    
    fft_conve = torchaudio.transforms.FFTConvolve()
    p_out     = fft_conve(sig, HRIR)
    
    return p_out

def plot_pyfar(Hnm,Ynm,fs,sig_path,save_path,ang,bp):
    dry_sig = pf.io.read_audio(sig_path)
    HRIR   = f.spatial_interpolation(Hnm,Ynm,True).cpu().numpy()
    HRIR = np.transpose(HRIR,[0,2,1])
    print("HRIR on lateral plane shape: ", HRIR.shape)
    HRIR_pf = pf.Signal(HRIR, fs)
    HRIR_pf_bp = pf.dsp.filter.butterworth(HRIR_pf, 8, bp, btype='bandpass')
    #HRIR_pf = pf.dsp.minimum_phase(HRIR_pf)
    
    p_out       = pf.dsp.convolve(HRIR_pf[ang],dry_sig)
    p_out       = pf.dsp.normalize(p_out,channel_handling='max')
    p_out_bp    = pf.dsp.convolve(HRIR_pf_bp[ang],dry_sig)
    p_out_bp    = pf.dsp.normalize(p_out_bp,channel_handling='max')
    print("The channel-wise energy: ", pf.dsp.energy(HRIR_pf[ang]))

    ax = pf.plot.time_freq(HRIR_pf[ang], label='120 deg HRIR')
    ax[0].legend()
    plt.show()
    
    ax = pf.plot.time_freq(HRIR_pf_bp[ang], label='HRIR band-pass')
    ax[0].legend()
    plt.show()

    save_name = save_path +f"p_{ang}.wav"
    dry_sig = pf.io.write_audio(p_out,save_name)
    save_name = save_path +f"p_BP_{ang}.wav"
    dry_sig = pf.io.write_audio(p_out_bp,save_name)
    

def save_as_sofa(Hnm,Ynm,SourcePosition,fs,save_path,HRIR_name):
    # Save the SH coeffitients as a sofa file (in [time x space x left/right])
    HRIR_IR = f.spatial_interpolation(Hnm,Ynm,True).cpu().numpy()
    HRIR_IR = np.transpose(HRIR_IR,[0,2,1])
    sofa = sf.Sofa("SimpleFreeFieldHRIR")
    sofa.Data_SamplingRate = fs
    sofa.SourcePosition = SourcePosition
    sofa.Data_IR = HRIR_IR
    sofa.delete("SourceUp")
    sf.write_sofa(save_path+HRIR_name, sofa)
    #data_ir, source_coordinates, receiver_coordinates = pf.io.read_sofa(save_path+HRIR_name)
    #index, *_ = source_coordinates.find_nearest_k(90, 0, 3.25, k=1, domain='sph', convention='top_elev', unit='deg', show=True)
    #_, mask = source_coordinates.find_slice('elevation', unit='deg', value=0, show=True)



    
def save_audio_to_path(save_path,x,fs):
    # Iterate over the first dimension of tensor x
    for i in range(x.shape[0]):
        # Get the audio data for the current slice
        audio_data = x[i].cpu()  # Assuming you're working on CPU
        
        # Construct the file path for saving the audio file
        file_path = save_path + f"audio_{i}.wav"
        
        # Write the audio data to the WAV file
        torchaudio.save(file_path, audio_data, fs, channels_first=True)
    
        #print(f"Audio file saved: {file_path}")
    
def import_matlab_data(data_path,device,shutup):
    # Import the Variables See Discriptions above for more info
    #print(os.path.basename(data_path))
    # Assign the local names for the variables (extract theme from the matlab dict)
    mat     = scipy.io.loadmat(data_path)
    fs      = mat["fs"][0,0]
    f_band  = mat["f_band"][0]
    N_low   = mat["N_low"][0]
    N_high  = mat["N_high"][0]
    f_vec   = mat["f_vec"]
    ang_vec = mat["ang_vec"]
    omega      = mat["omega"]
    omega_az   = mat["omega_az"]
    nfft    = (f_vec.shape[0]-1)*2
    
    Hnm_high      = mat["Hnm_high"]
    Hnm_high      = torch.from_numpy(Hnm_high).to(device)
    Hnm_low       = mat["Hnm_low"]
    Hnm_low       = torch.from_numpy(Hnm_low).to(device)
    Hnm_mls       = mat["Hnm_mls"]
    Hnm_mls       = torch.from_numpy(Hnm_mls).to(device)
    
    Y_high_lebedev = mat["Y_high_lebedev"]
    Y_high_lebedev = torch.from_numpy(Y_high_lebedev).to(device)
    Y_high_az      = mat["Y_high_az"]
    Y_high_az      = torch.from_numpy(Y_high_az).to(device)
    Y_low_lebedev  = mat["Y_low_lebedev"]
    Y_low_lebedev  = torch.from_numpy(Y_low_lebedev).to(device)
    Y_low_az       = mat["Y_low_az"]
    Y_low_az       = torch.from_numpy(Y_low_az).to(device)
    
    
    p_f_high_lebedev = mat["p_f_high_lebedev"]
    p_f_high_lebedev = torch.from_numpy(p_f_high_lebedev).to(device)
    ILD_ref          = mat["ILD_ref"]
    ILD_ref          = torch.from_numpy(ILD_ref).to(device)

    # This is a One time operation for getting the constanct that are needed for the AK-ILD calculations
    p_high_az_t  = f.spatial_interpolation(Hnm_high,Y_high_az,True)
    P            = torch.permute(p_high_az_t, (0, 2, 1))
    AK_Nmax,AK_f_c,AK_n,AK_C     = ak.AKerbILD_short_p1(P, f_band, fs)
    AK_C  = torch.from_numpy(AK_C).to(device)
    
    if not(shutup):
        print("\n\n----------Variables------------")
        print("p_f_high_lebedev",p_f_high_lebedev.shape, "\t",p_f_high_lebedev.dtype)
        print("ILD_ref \t",ILD_ref.shape, "\t\t",ILD_ref.dtype)
        print("")
        print("Hnm_high\t",Hnm_high.shape, "\t",Hnm_high.dtype)
        print("Hnm_low\t\t",Hnm_low.shape, "\t",Hnm_low.dtype)
        print("Hnm_mls\t\t",Hnm_mls.shape, "\t",Hnm_mls.dtype)
        print("")
        print("Y_high_lebedev\t",Y_high_lebedev.shape, "\t",Y_high_lebedev.dtype)
        print("Y_high_az\t",Y_high_az.shape, "\t",Y_high_az.dtype)
        print("Y_low_lebedev\t",Y_low_lebedev.shape, "\t\t",Y_low_lebedev.dtype)
        print("Y_low_az\t",Y_low_az.shape, "\t\t",Y_low_az.dtype)
        print("")
        
        
        
        print("f_band\t\t",f_band.shape, "\t\t\t\t",f_band.dtype)
        print("fs\t\t",fs.shape, "\t\t\t\t",fs.dtype)
        print("f_vec\t\t",f_vec.shape, "\t\t\t",f_vec.dtype)
        print("ang_vec\t\t",ang_vec.shape, "\t\t\t",ang_vec.dtype)
        print("N_high\t\t",N_high.shape, "\t\t\t\t",N_high.dtype)
        print("N_low\t\t",N_low.shape, "\t\t\t\t",N_low.dtype)

    p_ref_az_t = f.spatial_interpolation(Hnm_high,Y_high_az,True)

    return omega, omega_az, fs, f_band, N_low, N_high, f_vec, ang_vec, nfft, Hnm_high, Hnm_low, Hnm_mls, Y_high_lebedev, Y_high_az, Y_low_lebedev, Y_low_az, p_f_high_lebedev, ILD_ref, p_high_az_t, AK_Nmax, AK_f_c, AK_n, AK_C, p_ref_az_t
    



def plot_sofa_group_delay(path_sofa,figures_savepath,fig_name,is_save,shutup):
    ref, sources, *_ = pf.io.read_sofa(path_sofa)
    horizontal = np.abs(sources.elevation) < (.1 / 180 * np.pi)
    f_vec = np.linspace(0,ref.sampling_rate/2,ref.n_bins)
    impulse_response = torch.tensor(ref[horizontal, 0].time).permute(1,0) # take the left ear group delay
    group_delay_result_L = itd_est.group_delay_multi_channel(impulse_response)

    
    plt.figure
    plt.semilogx(f_vec,group_delay_result_L)
    plt.grid(True)
    plt.ylim(0,600)
    if is_save:
        plt.savefig(figures_savepath+ "/"+fig_name)
    if not(shutup):
        plt.show()
    else:
        plt.close()  


# Function to print a loading bar
def print_loading_bar(progress,ang):
    bar_length = 29
    ang_az = ang[0]
    ang_el = ang[1]
    filled_length = int(progress * bar_length)
    bar = '=' * filled_length + '-' * (bar_length - filled_length)
    print(f'\r[{bar}] ({ang_az:03d},{ang_el:03d})[az,el]', end='', flush=True)


def gen_source_alurazation(angles,az_el_flag,Hnm,Ynm,Omega,dry_sig_path,save_path,save_name,fs,shutup):
    bp = [20, 20e3]
    dry_sig    = pf.io.read_audio(dry_sig_path)
    HRIR       = f.spatial_interpolation(Hnm,Ynm,True).cpu().numpy()
    HRIR       = np.transpose(HRIR,[0,2,1])
    HRIR_pf    = pf.Signal(HRIR, fs)
    HRIR_pf    = pf.dsp.filter.butterworth(HRIR_pf, 8, bp, btype='bandpass')

    # Convert Omega to pf.Coordinates
    coardinate = pf.Coordinates.from_spherical_colatitude(Omega[:,1],Omega[:,0],Omega[:,2]) # [az, el, r]

    # Get the the idecies of the lateral of median plane (azimuth == 0 or 180)    
    if az_el_flag:
        mask = np.squeeze(np.logical_or(coardinate.azimuth == 0, coardinate.azimuth == np.pi))
        coardinate_new_plane = coardinate[mask]
        sort = np.argsort(coardinate_new_plane.polar)
        coardinate_new_plane = coardinate_new_plane[sort]
        # Convert angles to pf.Coordinates
        angles = np.deg2rad(angles)
        az_column = np.full((angles.shape[0], 1),0)
        radius_column = np.full((angles.shape[0], 1), Omega[0,2])
        # Stack the original matrix with the constant column horizontally
        angles = np.hstack((az_column, angles))
        angles = np.hstack((angles, radius_column))
        to_find = pf.Coordinates.from_spherical_side(angles[:,0], angles[:,1],  angles[:,2])
        # find the first and last indices in coardinate_new_plane acodring to angles
        index, distance = coardinate_new_plane.find_nearest(to_find)
        
    else:
        tol = np.deg2rad(0.5)
        mask = np.squeeze(np.logical_and(coardinate.elevation <= (0+tol),coardinate.elevation >= (0-tol)))
        coardinate_new_plane = coardinate[mask]
        sort = np.argsort(coardinate_new_plane.azimuth)
        coardinate_new_plane = coardinate_new_plane[sort]
        # Convert angles to pf.Coordinates
        angles = np.deg2rad(angles)
        el_column = np.full((angles.shape[0], 1),np.pi/2)
        radius_column = np.full((angles.shape[0], 1), Omega[0,2])
        # Stack the original matrix with the constant column horizontally
        angles = np.hstack((angles, el_column))
        angles = np.hstack((angles, radius_column))
        to_find = pf.Coordinates.from_spherical_colatitude(angles[:,0], angles[:,1],  angles[:,2])
        # find the first and last indices in coardinate_new_plane acodring to angles
        index, distance = coardinate_new_plane.find_nearest(to_find)
    

    HRIR_pf = HRIR_pf[mask,:]
    HRIR_pf = HRIR_pf[sort]

    
    
    if not(shutup):
        coardinate_new_plane.show(index)
        plt.show()
    
    num_of_HRIRs = index[0][1] - index[0][0] +1
  

    # using 50% overlap between two adjacent points using triangular window ###
   
    stride = np.ceil(dry_sig.n_samples / (num_of_HRIRs+1)).astype(int)
    window_len_samples = 2*stride + 1

    add_zeros = stride*(num_of_HRIRs+1) - dry_sig.n_samples
    # Zero pad the dry signal so the windows will fit
    padded = np.zeros((1,1+stride*(num_of_HRIRs+1)))
    padded[:, :dry_sig.time.shape[1]] = dry_sig.time
    dry_sig = pf.Signal(padded, dry_sig.sampling_rate)

    # initialize output signal
    output_signal = pf.Signal(np.zeros(dry_sig.n_samples+HRIR_pf.n_samples), dry_sig.sampling_rate)

    
    # create a triangular window
    window = np.bartlett(window_len_samples)

    # the length of the signals in term of window length
    signal_frames_len = int(output_signal.n_samples / window_len_samples)


    # initialize output signal
    output_signal = pf.Signal(np.zeros(dry_sig.n_samples+HRIR_pf.n_samples), dry_sig.sampling_rate)


    
    for frame_idx in range(num_of_HRIRs):
        first_sample = frame_idx*stride
        last_sample  = frame_idx*stride+window_len_samples
        curr_corr_index = index[0][0] + frame_idx
        curr_ang = [np.round(np.rad2deg(coardinate_new_plane.azimuth[curr_corr_index])).astype(int),np.round(np.rad2deg(coardinate_new_plane.elevation[curr_corr_index])).astype(int)]


        
        curr_window = np.pad(window, (first_sample, dry_sig.n_samples - last_sample), 'constant', constant_values=(0, 0))
        if frame_idx == 0:
            curr_window[:stride] = 1
        elif frame_idx == num_of_HRIRs - 1:
            curr_window[-stride:] = 1
            
        curr_window_pf = pf.Signal(curr_window, dry_sig.sampling_rate)
        dry_sig_window = pf.multiply((dry_sig, curr_window_pf), domain='time')

        '''
        plt.figure()
        pf.plot.time(curr_window_pf)
        plt.show()
        '''

        if frame_idx == 0:
            if not(shutup):
                print("start angle: ",curr_ang," [deg]\n")
            output_signal = pf.dsp.convolve(HRIR_pf[curr_corr_index],dry_sig_window)
        else:
            output_signal = output_signal + pf.dsp.convolve(HRIR_pf[curr_corr_index],dry_sig_window)
        # Simulate progress
        print_loading_bar(((frame_idx+1) / num_of_HRIRs),curr_ang)
        #print_loading_bar((frame_idx / num_of_HRIRs),4)
        time.sleep(0.01)  # Simulate some work being done


    if not(shutup):
        print("\n\nend angle: ",curr_ang,"[deg]")
    sig_length_time = output_signal.n_samples / fs
    movmentspeed = (num_of_HRIRs) / sig_length_time
    movmentspeed_in_sample = (num_of_HRIRs) /  output_signal.n_samples



    if not(shutup):
        print("signal time: ", sig_length_time," seconds")
        print("angle speed: ", movmentspeed," deg per sec")
        print("total frame number: ", num_of_HRIRs," (including 50% overlap)")
        print("num of sampels in each frame: ", window_len_samples," [samples]")
        print("\n Saving to disk")
    
    output_signal    = pf.dsp.normalize(output_signal,channel_handling='max')
    
    tmp = save_path+"/"+save_name+".wav"
    pf.io.write_audio(output_signal,tmp)
    if not(shutup):
        print("\n Done...\n\n")


def get_static_source(ang,Hnm,Ynm,Omega,dry_sig_path,save_path,save_name,fs,shutup):
    bp = [20, 20e3]
    dry_sig    = pf.io.read_audio(dry_sig_path)
    dry_sig    = pf.dsp.filter.butterworth(dry_sig, 8, bp, btype='bandpass')
    HRIR       = f.spatial_interpolation(Hnm,Ynm,True).cpu().numpy()
    HRIR       = np.transpose(HRIR,[0,2,1])
    HRIR_pf    = pf.Signal(HRIR, fs)
    #HRIR_pf    = pf.dsp.filter.butterworth(HRIR_pf, 8, bp, btype='bandpass')

    # Convert Omega to pf.Coordinates
    coardinate = pf.Coordinates.from_spherical_colatitude(Omega[:,1],Omega[:,0],Omega[:,2]) # [az, el, r]
    ang = np.deg2rad(ang)
    radius_column = np.full((ang.shape[0], 1), Omega[0,2])
    # Stack the original matrix with the constant column horizontally
    ang = np.hstack((ang, radius_column))
    to_find = pf.Coordinates.from_spherical_colatitude(ang[:,0], ang[:,1],  ang[:,2])
    # find the index in coardinate to angles
    index, distance = coardinate.find_nearest(to_find)
    if not(shutup):
        coardinate.show(index)
        plt.show()
    output_signal = pf.dsp.convolve(HRIR_pf[index],dry_sig)
    if not(shutup):
        print("\n Saving to disk")
    output_signal    = pf.dsp.normalize(output_signal,channel_handling='max')
    tmp = save_path+"/"+save_name+".wav"
    pf.io.write_audio(output_signal,tmp)
    if not(shutup):
        print("\n Done...\n\n")
    
    
    
def get_sound_aluriazitions(Hnm,Y,omega,ang,results_save_path,dry_sig_name,base_dir,save_name,fs,shutup):
    #dry_sig_name = "casta.wav"
    radius = 3.25
    # Create the third column with a constant value of 3.25
    radius = np.full((omega.shape[0],1), radius)
    # Stack the original matrix with the constant column horizontally
    new_omega = np.hstack((omega, radius))
    
    save_data_path = os.path.join(results_save_path, "Auralization")
    save_data_path =  os.path.join(save_data_path, dry_sig_name.split(".")[0])
    dry_sig_path = os.path.join(base_dir, "dry_signals", dry_sig_name)
    median_flag = False # the median of lateral plane flag, if True then median of false then lateral
    
    if median_flag:
        save_data_path_aura =  os.path.join(save_data_path, "Median_Plane")
    else:
        save_data_path_aura =  os.path.join(save_data_path, "Horizontal_Plane")
        
    srat_stop_angles = np.atleast_2d(np.array([0,300])).T; # if lateral then 0:360 where 90 is the left ear, if madian the use polar -90:270 where 0 is front and -90 is down, 90 is up and 180 is back
    
    if not os.path.isdir(save_data_path_aura):
            os.makedirs(save_data_path_aura)
    
   
    gen_source_alurazation(srat_stop_angles,median_flag,Hnm,Y,new_omega,dry_sig_path,save_data_path_aura,save_name,fs,shutup)
    
    
    median_flag = True # the median of lateral plane flag, if True then median of false then lateral
    
    if median_flag:
        save_data_path_aura =  os.path.join(save_data_path, "Median_Plane")
    else:
        save_data_path_aura =  os.path.join(save_data_path, "Lateral_Plane")
    srat_stop_angles = np.atleast_2d(np.array([-80,180])).T; # if lateral then 0:360 where 90 is the left ear, if madian the use polar -90:270 where 0 is front and -90 is down, 90 is up and 180 is back

    
    if not os.path.isdir(save_data_path_aura):
            os.makedirs(save_data_path_aura)
    
    gen_source_alurazation(srat_stop_angles,median_flag,Hnm,Y,new_omega,dry_sig_path,save_data_path_aura,save_name,fs,shutup)




    foo = "Static_" + np.array2string(ang[0])
    save_data_path_static =  os.path.join(save_data_path, foo)
    
    
    
    if not os.path.isdir(save_data_path_static):
            os.makedirs(save_data_path_static)
    
    
  
    get_static_source(ang,Hnm,Y,new_omega,dry_sig_path,save_data_path_static,save_name,fs,shutup)
    

def plot_time_freq_example(path_sofa,ang,fs,figures_savepath,file_name,is_save,shutup):
    ref, sources, *_ = pf.io.read_sofa(path_sofa)

    ang = np.deg2rad(ang)
    radius = 3.25
    radius_column = np.full((ang.shape[0], 1), radius)
    # Stack the original matrix with the constant column horizontally
    ang = np.hstack((ang, radius_column))

    to_find = pf.Coordinates.from_spherical_colatitude(ang[:,0], ang[:,1],  ang[:,2])
    # find the index in coardinate to angles
    index, distance = sources.find_nearest(to_find)


    HRIR_pf = ref[index]
    

    ax = pf.plot.time_freq(HRIR_pf)
    #ax[0].legend()
    if is_save:
        plt.savefig(figures_savepath+ file_name)
    if not(shutup):
        plt.show()
    else:
        plt.close()
 



    

def start(data_path,results_save_path,epochs=300,lambda_vec=[1,1,1,1],dry_sig_name="casta.wav",flag_s2=False,shutup=True,is_save=False):

    figures_savepath = os.path.join( results_save_path,"training_figures")

    # Setup device
    has_gpu = torch.cuda.is_available()
    has_mps = torch.backends.mps.is_built()
    device = "mps" if torch.backends.mps.is_built() \
        else "gpu" if torch.cuda.is_available() else "cpu"
    if not(shutup):
        print(f"Python Platform: {platform.platform()}")
        print(f"PyTorch Version: {torch.__version__}")
        print()
        print(f"Python {sys.version}")
        print("NVIDIA/CUDA GPU is", "available" if has_gpu else "NOT AVAILABLE")
        print("MPS (Apple Metal) is", "AVAILABLE" if has_mps else "NOT AVAILABLE")
    device = "cpu"
    if not(shutup):
        print(f"Target device is {device}")

    if is_save:
        if not os.path.isdir(figures_savepath):
            os.makedirs(figures_savepath)

    omega, omega_az, fs, f_band, N_low, N_high, f_vec, ang_vec, nfft, Hnm_high, Hnm_low, Hnm_mls, Y_high_lebedev, Y_high_az, Y_low_lebedev, Y_low_az, p_f_high_lebedev, ILD_ref, p_high_az_t, AK_Nmax, AK_f_c, AK_n, AK_C, p_ref_az_t = import_matlab_data(data_path,device,shutup)


    # Source position for the sofa export
    SourcePosition  = np.rad2deg(omega)
    source_distance = 3.25
    source_distance_column = np.full((SourcePosition.shape[0], 1), source_distance)
    SourcePosition = np.concatenate((SourcePosition, source_distance_column), axis=1)
    SourcePosition[:,0] = -1*(SourcePosition[:,0] - 90)  # change from [0,180) to [-90,90) with +30 left and -30 right 0 is facing forword
    SourcePosition[:, [0, 1]] = SourcePosition[:, [1, 0]]

    save_path_sofa = os.path.join(results_save_path, "sofa_export")
    if not os.path.isdir(save_path_sofa):
            os.makedirs(save_path_sofa)

    # Save the reference sofa and the MagLS sofa for the baumgartner2014 model
    if is_save:
        HRIR_name = "/ref.sofa"
        save_as_sofa(Hnm_high,Y_high_lebedev,SourcePosition,int(fs),save_path_sofa,HRIR_name)
        
        path_sofa = save_path_sofa + HRIR_name
        fig_name = "left_ear_group_delay-reference.png"
        plot_sofa_group_delay(path_sofa,figures_savepath,fig_name,is_save,shutup)
        
        HRIR_name = "/MagLS.sofa"
        save_as_sofa(Hnm_mls,Y_low_lebedev,SourcePosition,int(fs),save_path_sofa,HRIR_name)
        path_sofa = save_path_sofa + HRIR_name
        fig_name = "left_ear_group_delay-MagLS.png"
        plot_sofa_group_delay(path_sofa,figures_savepath,fig_name,is_save,shutup)
    
    

    model    = NN_v2_simp(Hnm_low.shape[0],Hnm_low.shape[1])
    model.to(device)
    
    optimizer = torch.optim.SGD(model.parameters(), lr=0.00003, momentum=0.9)
    optimizer.zero_grad()

    x = Hnm_mls  # initial solution is Ambisonics MagLS
    #x = Hnm_low # initial solution is Ambisonics LS
    x.to(device)

    if not(shutup):  
        print("\n\n-----------NN summary---------")
        result = summary(model,input_size= x.shape,dtypes=[torch.complex128],verbose=2,col_width=13,
            col_names=["kernel_size","input_size", "output_size", "num_params", "mult_adds"], row_settings=["var_names"],)

    h_nmse = []
    h_mag = []
    h_ILD = []
    h_total = []
    h_diff = []
    h_color = []
    h_ITD = []
    h_qe = []
    h_pe = []
    h_color_weighted = []
    h_mag_weighted = []


    

    cutoff_freq = 3e3
    #cutoff_freq = 1.5e3
    trashhold   = -20
    x_lr = f.spatial_interpolation(Hnm_high,Y_high_az,True).float()
    ITD_ref = itd_est.get_ITD_group_delay(x_lr,fs,cutoff_freq,False) #returns the micro sec ITD usinng the group delay method
    GD_ref   = itd_est.get_muilti_ch_group_delay(x_lr,fs, cutoff_freq,False) # returns the group delay in [freq x space x ears]


    
    if flag_s2:
        HRIR_ref = f.spatial_interpolation(Hnm_high,Y_high_lebedev,True)
        ILD_ref_opt = ak.clc_ILD(HRIR_ref, AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)
    else:
        HRIR_ref = f.spatial_interpolation(Hnm_high,Y_high_az,True)
        ILD_ref_opt = ak.clc_ILD(HRIR_ref, AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft)
    
    pos_grad_ref_fc = calc_positive_grad(Hnm_high,Y_high_lebedev,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft,True)

    pos_grad_ref_f  = calc_positive_grad_pyfar(Hnm_high,Y_high_lebedev,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft,True)


    

    #pos_grad_ref_bins = calc_positive_grad_bins(Hnm_high,Y_high_az,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft,False)

    #pos_grad_ref_bins = calc_positive_grad_interpulate(Hnm_high,Y_high_az,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C, nfft,False)
    

    for epoch in tqdm (range (epochs), desc="Loading..."):

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        output = model(x)
        output_dict   = loss_function_v4(output,Y_low_lebedev,Y_low_az,p_f_high_lebedev,p_ref_az_t,
                                         ILD_ref_opt,ITD_ref,GD_ref,pos_grad_ref_fc,pos_grad_ref_f,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C,nfft,
                                         lambda_vec,save_path_sofa,SourcePosition,flag_s2,norm_flag=False)

        output_dict["e_total"].backward(retain_graph=True)
        optimizer.step()
        
        
        h_color_weighted     = np.append(h_color_weighted,output_dict["e_color_weighted"].detach().cpu().numpy()*lambda_vec[6])
        h_mag_weighted       = np.append(h_mag_weighted,torch.norm(output_dict["e_mag_weighted"],dim=0).detach().cpu().numpy()*lambda_vec[7])
        
        h_qe       = np.append(h_qe,output_dict["qe"])
        h_pe       = np.append(h_pe,output_dict["pe"])
        
        h_ITD       = np.append(h_ITD,output_dict["e_ITD"].detach().cpu().numpy()*lambda_vec[5])
        h_color     = np.append(h_color,output_dict["e_color"].detach().cpu().numpy()*lambda_vec[4])
        h_diff      = np.append(h_diff,output_dict["e_mag_diff"].detach().cpu().numpy()*lambda_vec[3])
        
        h_ILD       = np.append(h_ILD,torch.norm(output_dict["e_ILD"],dim=0).detach().cpu().numpy()*lambda_vec[0])
        
        h_nmse      = np.append(h_nmse,torch.norm(output_dict["e_nmse"],dim=0).detach().cpu().numpy()*lambda_vec[1])
        h_mag       = np.append(h_mag,torch.norm(output_dict["e_mag"],dim=0).detach().cpu().numpy()*lambda_vec[2])
        h_total     = np.append(h_total,torch.norm(output_dict["e_total"],dim=0).detach().cpu().numpy())
        

    
        #print(h_color)
        #print(h_diff)
        #print(h_ILD)
        #print(h_nmse)
        #print(h_mag)
        

    if not(shutup):    
        print('Finished Training')

    
    """
    max_value = np.max(h_nmse, axis=0)
    normalized_h_nmse = h_nmse / max_value
    
    max_value = np.max(h_mag, axis=0)
    normalized_h_mag = h_mag / max_value
    
    max_value = np.max(h_ILD, axis=0)
    normalized_h_ILD = h_ILD / max_value
    
    max_value = np.max(h_diff, axis=0)
    normalized_h_diff = h_diff / max_value

    max_value = np.max(h_color, axis=0)
    normalized_h_color = h_color / max_value
    
    max_value = np.max(h_ITD, axis=0)
    normalized_h_ITD = h_ITD / max_value

    max_value = np.max(h_color_weighted, axis=0)
    normalized_h_color_weighted = h_color_weighted / max_value
    
    max_value = np.max(h_mag_weighted, axis=0)
    normalized_h_mag_weighted = h_mag_weighted / max_value
    """
    normalized_h_nmse  = h_nmse
    normalized_h_mag   = h_mag
    normalized_h_ILD   = h_ILD
    normalized_h_diff  = h_diff
    normalized_h_color = h_color
    normalized_h_ITD   = h_ITD
    normalized_h_color_weighted = h_color_weighted
    normalized_h_mag_weighted = h_mag_weighted
    
    plt.plot(normalized_h_nmse)
    plt.plot(normalized_h_mag)
    plt.plot(normalized_h_ILD)
    plt.plot(normalized_h_diff)
    plt.plot(normalized_h_color)
    plt.plot(normalized_h_ITD)
    plt.plot(normalized_h_color_weighted)
    plt.plot(normalized_h_mag_weighted)
    plt.grid()
    plt.xlabel("iterations")
    plt.ylabel("error")
    plt.title("Training Curves")
    plt.legend(["e NMSE","e Magnitude","e ILD","e_mag_diff","e_color","e_ITD","e_color_weighted","e_mag_weighted"]);
    
    if is_save:
        plt.savefig(figures_savepath+ "/Training_Curves.png")
    
    if not(shutup):
        plt.show()
    else:
        plt.close() 
    
    

    # Save the initial and final values 
    # Define the arrays
    h1 = [h_nmse[0], h_nmse[-1]]
    h2 = [h_mag[0], h_mag[-1]]
    h3 = [h_ILD[0], h_ILD[-1]]
    h4 = [h_diff[0], h_diff[-1]]
    h5 = [h_color[0], h_color[-1]]
    h6 = [h_ITD[0], h_ITD[-1]]
    h7 = [h_color_weighted[0], h_color_weighted[-1]]
    h8 = [h_mag_weighted[0], h_mag_weighted[-1]]
    
    # Define the file path
    file_path = figures_savepath + "/error_values.txt"
    
    # Define array names and values
    arrays = [("nmse error [dB]", h1), ("magnitude error [dB]", h2), ("ILD error [dB]", h3),("diff error", h4),("color error", h5),("ITD error [samples]", h6),("color weighted error ", h7),("magnitude weighted error [dB]", h8)]
    
    # Write the arrays to the file
    with open(file_path, "w") as file:
        for name, values in arrays:
            file.write(f"{name}:\n")
            file.write(f"{values[0]}, {values[1]}\n\n")
    
    plt.plot(h_qe)
    plt.plot(h_pe)
    plt.grid()
    plt.xlabel("iterations")
    plt.ylabel("error")
    plt.title("Training Curves")
    plt.legend(["Quadrant errors (%)","Local polar RMS error (deg)"]);
    if is_save:
        plt.savefig(figures_savepath+ "/qe_pe_curves.png")
    if not(shutup):
        plt.show()
    else:
        plt.close()

    Hnm_imls = output.detach()

    x_lr = f.spatial_interpolation(Hnm_imls,Y_low_az,True).float()
    #ITD_imls = itd_est.get_ITD_onset_threshold(x_lr,fs,trashhold,cutoff_freq,False) #returns the micro sec ITD usinng the onset method
    ITD_imls = itd_est.get_ITD_group_delay(x_lr,fs,cutoff_freq,False) #returns the micro sec ITD usinng the onset method
    x_lr = f.spatial_interpolation(Hnm_mls,Y_low_az,True).float()
    #ITD_mls = itd_est.get_ITD_onset_threshold(x_lr,fs,trashhold,cutoff_freq,False) #returns the micro sec ITD usinng the onset method
    ITD_mls = itd_est.get_ITD_group_delay(x_lr,fs,cutoff_freq,False) #returns the micro sec ITD usinng the onset method


    ITD_ref_tmp = ITD_ref.detach().numpy()
    ITD_ref_tmp = (ITD_ref_tmp/fs)*1e6
    ITD_imls = ITD_imls.detach().numpy()
    ITD_imls = (ITD_imls/fs)*1e6
    ITD_mls = ITD_mls.detach().numpy()
    ITD_mls = (ITD_mls/fs)*1e6
    
    plt.figure
    plt.plot(ITD_ref_tmp, label='Reference')
    plt.plot(ITD_mls, label='MagLS')
    plt.plot(ITD_imls, label='iMagLS')
    plt.xlabel('Incident angle')
    plt.ylabel('Value [micro sec]')
    plt.title('ITD')
    plt.grid(True)
    plt.legend()
    if is_save:
        plt.savefig(figures_savepath+ "/ITD_curves.png")
    if not(shutup):
        plt.show()
    else:
        plt.close()
        
    # Reassign the original low frequency coefficients
    idx_min = (np.abs(f_vec - 2e3)).argmin()
    #Hnm_imls[:,0:idx_min,:] = Hnm_mls[:,0:idx_min,:]

    
    # Calc some plots
    lambda_vec_tests = [1,1,1,1,1,1,0,0]
    res_imls   = loss_function_v4(Hnm_imls,Y_low_lebedev,Y_low_az,p_f_high_lebedev,p_ref_az_t,
                                     ILD_ref,ITD_ref,GD_ref,pos_grad_ref_fc,pos_grad_ref_f,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C,nfft,
                                     lambda_vec_tests,save_path_sofa,SourcePosition,flag_s2=False,norm_flag=True)
    
    res_ls   = loss_function_v4(Hnm_low,Y_low_lebedev,Y_low_az,p_f_high_lebedev,p_ref_az_t,
                                     ILD_ref,ITD_ref,GD_ref,pos_grad_ref_fc,pos_grad_ref_f,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C,nfft,
                                     lambda_vec_tests,save_path_sofa,SourcePosition,flag_s2=False,norm_flag=True)
    
    res_mls   = loss_function_v4(Hnm_mls,Y_low_lebedev,Y_low_az,p_f_high_lebedev,p_ref_az_t,
                                     ILD_ref,ITD_ref,GD_ref,pos_grad_ref_fc,pos_grad_ref_f,AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C,nfft,
                                     lambda_vec_tests,save_path_sofa,SourcePosition,flag_s2=False,norm_flag=True)
    

    e_ILD_imls  = res_imls["e_ILD"].detach().cpu().numpy()
    e_nmse_imls = res_imls["e_nmse"].detach().cpu().numpy()
    e_mag_imls  = res_imls["e_mag"].detach().cpu().numpy()
    ILD_imls    = res_imls["ILD"].detach()

    e_ILD_non_av_ls = torch.abs(ILD_ref - res_ls["ILD"])
    e_ILD_non_av_mls = torch.abs(ILD_ref - res_mls["ILD"])
    e_ILD_non_av_imls = torch.abs(ILD_ref - res_imls["ILD"])


    plt.plot(f_vec,res_ls["e_nmse"].detach().numpy())
    plt.plot(f_vec,res_mls["e_nmse"].detach().numpy())
    plt.plot(f_vec,res_imls["e_nmse"].detach().numpy())
    plt.plot(f_vec,res_ls["e_mag"].detach().numpy())
    plt.plot(f_vec,res_mls["e_mag"].detach().numpy())
    plt.plot(f_vec,res_imls["e_mag"].detach().numpy())
    plt.xlim([99,20e3])
    plt.xscale('log')
    plt.xlabel("f[Hz]")
    plt.grid(which='both')
    plt.title("Avaraged per ear NMSE and Magnitude Errors")
    plt.legend(["LS NMSE","MagLS NMSE","iMagLS NMSE","LS magnitude","MagLS magnitude","iMagLS magnitude"])
    plt.ylabel("error[dB]");
    if is_save:
        plt.savefig(figures_savepath+ "/freq_errors.png")
    if not(shutup):
        plt.show()
    else:
        plt.close()
    
    

    plt.plot(ang_vec,torch.mean(ILD_ref,dim=0))
    plt.plot(ang_vec,torch.mean(res_ls["ILD"],dim=0).detach().numpy())
    plt.plot(ang_vec,torch.mean(res_mls["ILD"],dim=0).detach().numpy())
    plt.plot(ang_vec,torch.mean(ILD_imls,dim=0))
    plt.grid()
    plt.xlabel("azimuth [deg]")
    plt.ylabel("ILD value [dB]")
    plt.legend(["Reference","LS","MagLS","iMagLS"])
    plt.xlim((ang_vec[0], ang_vec[-1]))
    plt.title("ILD Value average over all central frequencies");
    if is_save:
        plt.savefig(figures_savepath+ "/ILD_Curves.png")
    if not(shutup):
        plt.show()
    else:
        plt.close()
    
    plt.plot(ang_vec,res_ls["e_ILD"].detach().numpy())
    plt.plot(ang_vec,res_mls["e_ILD"].detach().numpy())
    plt.plot(ang_vec,e_ILD_imls)
    plt.grid()
    plt.xlabel("azimuth [deg]")
    plt.ylabel("Error [dB]")
    plt.legend(["LS","MagLS","iMagLS"])
    plt.xlim((ang_vec[0], ang_vec[-1]))
    plt.title("ILD error average over all central frequencies");
    if is_save:
        plt.savefig(figures_savepath+ "/ILD_errors_ang.png")
    if not(shutup):
        plt.show()
    else:
        plt.close()


    plt.plot(AK_f_c,torch.mean(e_ILD_non_av_ls,dim=1).detach().numpy(), '--o')
    plt.plot(AK_f_c,torch.mean(e_ILD_non_av_mls,dim=1).detach().numpy(), '--o')
    plt.plot(AK_f_c,torch.mean(e_ILD_non_av_imls,dim=1).detach().numpy(), '--o')
    plt.grid()
    plt.xlabel("f[Hz]")
    plt.ylabel("ILD error [dB]")
    plt.legend(["LS","MagLS","iMagLS"])
    plt.xlim((AK_f_c[0], AK_f_c[-1]))
    plt.title("ILD Error average over all lateral Directions");
    if is_save:
        plt.savefig(figures_savepath+ "/ILD_errors_freq.png")
    if not(shutup):
        plt.show()
    else:
        plt.close()


    
    e_color_fc_lr_space_imls = res_imls["e_color_fc_lr_space"].detach().cpu().numpy()
    e_color_fc_lr_space_mls  = res_mls["e_color_fc_lr_space"].detach().cpu().numpy()
    e_color_fc_lr_space_ls   = res_ls["e_color_fc_lr_space"].detach().cpu().numpy()


    e_color_fc_imls = np.squeeze(np.mean(e_color_fc_lr_space_imls,axis = 1))[:,0] # mean over space and take the left ear
    e_color_fc_mls = np.squeeze(np.mean(e_color_fc_lr_space_mls,axis = 1))[:,0] # mean over space and take the left ear
    e_color_fc_ls = np.squeeze(np.mean(e_color_fc_lr_space_ls,axis = 1))[:,0] # mean over space and take the left ear

    e_color_fc_imls_right = np.squeeze(np.mean(e_color_fc_lr_space_imls,axis = 1))[:,1] # mean over space and take the left ear
    e_color_fc_mls_right  = np.squeeze(np.mean(e_color_fc_lr_space_mls,axis = 1))[:,1] # mean over space and take the left ear
    e_color_fc_ls_right   = np.squeeze(np.mean(e_color_fc_lr_space_ls,axis = 1))[:,1] # mean over space and take the left ear


    plt.plot(AK_f_c,e_color_fc_ls, '--o')
    plt.plot(AK_f_c,e_color_fc_mls, '--o')
    plt.plot(AK_f_c,e_color_fc_imls, '--o')
    plt.plot(AK_f_c,e_color_fc_ls_right, '--o')
    plt.plot(AK_f_c,e_color_fc_mls_right, '--o')
    plt.plot(AK_f_c,e_color_fc_imls_right, '--o')
    plt.grid()
    plt.xlabel("f[Hz]")
    plt.ylabel("Coloration error [dB]")
    plt.legend(["LS-left","MagLS-left","iMagLS-left","LS-right","MagLS-right","iMagLS-right"])
    plt.xlim((AK_f_c[0], AK_f_c[-1]))
    plt.title("Coloration Error average over Space");
    if is_save:
        plt.savefig(figures_savepath+ "/color_errors_freq.png")
    if not(shutup):
        plt.show()
    else:
        plt.close()
    


    if is_save:
        HRIR_name = "/iMagLS.sofa"
        save_as_sofa(Hnm_imls,Y_low_lebedev,SourcePosition,int(fs),save_path_sofa,HRIR_name)
    

    path_sofa = save_path_sofa + HRIR_name
    fig_name = "left_ear_group_delay-iMagLS.png"
    plot_sofa_group_delay(path_sofa,figures_savepath,fig_name,is_save,shutup)


    # Save .wav files of results
    
    base_dir = os.path.join('..')
    save_name = "iMagLS"
    ang = np.atleast_2d(np.array([-60,90])) # (1,2) ang[0,0]=azimuth and ang[0,1]=elevation (spherical elevation system)

    get_sound_aluriazitions(Hnm_imls,Y_low_lebedev,omega,ang,results_save_path,dry_sig_name,base_dir,save_name,fs,shutup)

    figures_savepath_tf_plots = os.path.join(figures_savepath, "time_freq_plots")
    if not os.path.isdir(figures_savepath_tf_plots):
        os.makedirs(figures_savepath_tf_plots)

    
    ang = np.atleast_2d(np.array([-60,90])) # (1,2) ang[0,0]=azimuth and ang[0,1]=elevation (spherical elevation system)
    path_sofa = save_path_sofa + "/iMagLS.sofa"
    file_name = "/example_"  + np.array2string(ang[0]) + "_time_freq_plot-iMagLS.png"
    plot_time_freq_example(path_sofa,ang,fs,figures_savepath_tf_plots,file_name,is_save,shutup)
    path_sofa = save_path_sofa + "/MagLS.sofa"
    file_name = "/example_"  + np.array2string(ang[0]) + "_time_freq_plot-MagLS.png"
    plot_time_freq_example(path_sofa,ang,fs,figures_savepath_tf_plots,file_name,is_save,shutup)
    path_sofa = save_path_sofa + "/ref.sofa"
    file_name = "/example_"  + np.array2string(ang[0]) + "_time_freq_plot-ref.png"
    plot_time_freq_example(path_sofa,ang,fs,figures_savepath_tf_plots,file_name,is_save,shutup)

    ang = np.atleast_2d(np.array([0,90])) # (1,2) ang[0,0]=azimuth and ang[0,1]=elevation (spherical elevation system)
    path_sofa = save_path_sofa + "/iMagLS.sofa"
    file_name = "/example_"  + np.array2string(ang[0]) + "_time_freq_plot-iMagLS.png"
    plot_time_freq_example(path_sofa,ang,fs,figures_savepath_tf_plots,file_name,is_save,shutup)
    path_sofa = save_path_sofa + "/MagLS.sofa"
    file_name = "/example_"  + np.array2string(ang[0]) + "_time_freq_plot-MagLS.png"
    plot_time_freq_example(path_sofa,ang,fs,figures_savepath_tf_plots,file_name,is_save,shutup)
    path_sofa = save_path_sofa + "/ref.sofa"
    file_name = "/example_"  + np.array2string(ang[0]) + "_time_freq_plot-ref.png"
    plot_time_freq_example(path_sofa,ang,fs,figures_savepath_tf_plots,file_name,is_save,shutup)

    
                
    
    output_dict = {
        "Hnm_imls": Hnm_imls,
        "Hnm_mls": Hnm_mls,
        "Hnm_high": Hnm_high,
        "Y_high_az": Y_high_az,
        "Y_low_az": Y_low_az,
        "Y_high_lebedev": Y_high_lebedev,
        "Y_low_lebedev": Y_low_lebedev,
        "fs": fs,
        "omega": omega,
        "omega_az": omega_az
    }

    results_save_path_data = os.path.join(results_save_path, "data" )
    filename = "/data" + ".pkl"

    if is_save:
        if not os.path.isdir(results_save_path_data):
                os.makedirs(results_save_path_data)
        with open(results_save_path_data+filename, 'wb') as ff:
            pickle.dump(output_dict, ff)
            

    print("\n\n=======\nDone!\n=======")
    return output_dict