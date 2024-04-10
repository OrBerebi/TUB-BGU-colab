
import os
import sys
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
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import pyfar as pf
from torchinfo import summary




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
    e_mag = f.clc_e_mag(p_ref,p_lebedev,norm_flag)
    
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
    
    
def start(data_path,epochs=300,lambda_vec=[1,1,1,1],shutup=True,is_save=False):

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
    
    model    = NN_v2_simp(Hnm_low.shape[0],Hnm_low.shape[1])
    model.to(device)
    
    optimizer = torch.optim.SGD(model.parameters(), lr=0.00003, momentum=0.9)
    optimizer.zero_grad()

    x = Hnm_mls
    x.to(device)

    result = summary(model,input_size= x.shape,dtypes=[torch.complex128],verbose=2,col_width=13,
        col_names=["kernel_size","input_size", "output_size", "num_params", "mult_adds"], row_settings=["var_names"],)

    h_nmse = []
    h_mag = []
    h_ILD = []
    h_total = []
    h_enrgy = []


    for epoch in tqdm (range (epochs), desc="Loading..."):

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        output = model(x)
        loss, e_ILD_after, e_nmse_after, e_mag_after,ILD_after,e_enrgy_after   = loss_function_v2(output,Y_low_lebedev,
                                                                                    Y_low_az,p_f_high_lebedev,p_ref_az_t,ILD_ref,
                                                                                    AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C,
                                                                                    nfft,lambda_vec,norm_flag=False)

        loss.backward(retain_graph=True)
        optimizer.step()

        h_enrgy     = np.append(h_enrgy,e_enrgy_after.detach().cpu().numpy()*lambda_vec[3])
        e_ILD_after = torch.mean(torch.abs(ILD_ref - e_ILD_after),dim=0) #avarage over frequencies bands
        h_ILD       = np.append(h_ILD,torch.norm(e_ILD_after,dim=0).detach().cpu().numpy()*lambda_vec[0])
        h_nmse      = np.append(h_nmse,torch.norm(e_nmse_after,dim=0).detach().cpu().numpy()*lambda_vec[1])
        h_mag       = np.append(h_mag,torch.norm(e_mag_after,dim=0).detach().cpu().numpy()*lambda_vec[2])
        h_total     = np.append(h_total,torch.norm(loss,dim=0).detach().cpu().numpy())

    if not(shutup):    
        print('Finished Training')

    if not(shutup):
        plt.plot(h_nmse)
        plt.plot(h_mag)
        plt.plot(h_ILD)
        plt.plot(h_enrgy)
        plt.grid()
        plt.xlabel("iterations")
        plt.ylabel("error")
        plt.title("Training Curves")
        plt.legend(["e NMSE","e Magnitude","e ILD","e_mag_diff"]);
        if is_save:
            plt.savefig(data_path_dir+ "Training_Curves.png")
        plt.show()

    Hnm_imls = output.detach()

    # Reassign the original low frequency coefficients
    idx_min = (np.abs(f_vec - 2e3)).argmin()
    Hnm_imls[:,0:idx_min,:] = Hnm_mls[:,0:idx_min,:]

    
    # Calc some plots
    e_total_imls,e_ILD_imls, e_nmse_imls, e_mag_imls, ILD_imls = loss_function_v1(Hnm_imls,Y_low_lebedev,Y_low_az,p_f_high_lebedev,ILD_ref,
                                                                                     AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C,
                                                                                     nfft,lambda_vec,norm_flag=True)
    e_total_ls,e_ILD_ls, e_nmse_ls, e_mag_ls, ILD_ls = loss_function_v1(Hnm_low,Y_low_lebedev,Y_low_az,p_f_high_lebedev,ILD_ref,
                                                                                     AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C,
                                                                                     nfft,lambda_vec,norm_flag=True)
    e_total_mls,e_ILD_mls, e_nmse_mls, e_mag_mls, ILD_mls = loss_function_v1(Hnm_mls,Y_low_lebedev,Y_low_az,p_f_high_lebedev,ILD_ref,
                                                                                     AK_f_c, fs, f_band, AK_Nmax, AK_n, AK_C,
                                                                                     nfft,lambda_vec,norm_flag=True)



    e_ILD_imls  = e_ILD_imls.detach().cpu().numpy()
    e_nmse_imls = e_nmse_imls.detach().cpu().numpy()
    e_mag_imls  = e_mag_imls.detach().cpu().numpy()
    ILD_imls    = ILD_imls.detach()

    e_ILD_non_av_ls = torch.abs(ILD_ref - ILD_ls)
    e_ILD_non_av_mls = torch.abs(ILD_ref - ILD_mls)
    e_ILD_non_av_imls = torch.abs(ILD_ref - ILD_imls)


    
    if not(shutup):
        plt.plot(f_vec,e_nmse_ls)
        plt.plot(f_vec,e_nmse_mls)
        plt.plot(f_vec,e_nmse_imls)
        plt.plot(f_vec,e_mag_ls)
        plt.plot(f_vec,e_mag_mls)
        plt.plot(f_vec,e_mag_imls)
        plt.xlim([99,20e3])
        plt.xscale('log')
        plt.xlabel("f[Hz]")
        plt.grid(which='both')
        plt.title("Avaraged per ear NMSE and Magnitude Errors")
        plt.legend(["LS NMSE","MagLS NMSE","iMagLS NMSE","LS magnitude","MagLS magnitude","iMagLS magnitude"])
        plt.ylabel("error[dB]");
        plt.show()
        
        
    
        plt.plot(ang_vec,torch.mean(ILD_ref,dim=0))
        plt.plot(ang_vec,torch.mean(ILD_ls,dim=0))
        plt.plot(ang_vec,torch.mean(ILD_mls,dim=0))
        plt.plot(ang_vec,torch.mean(ILD_imls,dim=0))
        plt.grid()
        plt.xlabel("azimuth [deg]")
        plt.ylabel("ILD value [dB]")
        plt.legend(["Reference","LS","MagLS","iMagLS"])
        plt.xlim((ang_vec[0], ang_vec[-1]))
        plt.title("ILD Value average over all central frequencies");
        plt.show()
        
        plt.plot(ang_vec,e_ILD_ls)
        plt.plot(ang_vec,e_ILD_mls)
        plt.plot(ang_vec,e_ILD_imls)
        plt.grid()
        plt.xlabel("azimuth [deg]")
        plt.ylabel("Error [dB]")
        plt.legend(["LS","MagLS","iMagLS"])
        plt.xlim((ang_vec[0], ang_vec[-1]))
        plt.title("ILD error average over all central frequencies");
        plt.show()
    
    
        plt.plot(AK_f_c,torch.mean(e_ILD_non_av_ls,dim=1), '--o')
        plt.plot(AK_f_c,torch.mean(e_ILD_non_av_mls,dim=1), '--o')
        plt.plot(AK_f_c,torch.mean(e_ILD_non_av_imls,dim=1), '--o')
        plt.grid()
        plt.xlabel("f[Hz]")
        plt.ylabel("ILD error [dB]")
        plt.legend(["LS","MagLS","iMagLS"])
        plt.xlim((AK_f_c[0], AK_f_c[-1]))
        plt.title("ILD Error average over all lateral Directions");
        plt.show()

    
    


    if is_save:
        sig,fs_sig = torchaudio.load(sig_path)
        if not os.path.isdir(save_path+"/ref/"):
            os.makedirs(save_path+"/ref/")
        print("Binaural reproduction for reference")
        res = binaural_reproduction(Hnm_high,Y_high_az,sig)
        save_audio_to_path(save_path+"/ref/",res,fs_sig)
    
        if not os.path.isdir(save_path+"/LS/"):
            os.makedirs(save_path+"/LS/")
        print("Binaural reproduction for LS")
        res = binaural_reproduction(Hnm_low,Y_low_az,sig)
        save_audio_to_path(save_path+"/LS/",res,fs_sig)
    
        if not os.path.isdir(save_path+"/MagLS/"):
            os.makedirs(save_path+"/MagLS/")
        print("Binaural reproduction for MagLS")
        res = binaural_reproduction(Hnm_mls,Y_low_az,sig)
        save_audio_to_path(save_path+"/MagLS/",res,fs_sig)
    
        if not os.path.isdir(save_path+"/iMagLS/"):
            os.makedirs(save_path+"/iMagLS/")
        print("Binaural reproduction for iMagLS")
        res = binaural_reproduction(Hnm_imls,Y_low_az,sig)
        save_audio_to_path(save_path+"/iMagLS/",res,fs_sig)
    
    return Hnm_imls, Hnm_mls, Hnm_high, Y_high_az, Y_low_az, fs




    


    




    

