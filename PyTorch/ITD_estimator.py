import torch
import scipy.signal as signal
import matplotlib.pyplot as plt

def get_onset(x,onset_threshold_dB):
    x_abs = torch.abs(x)
    onset_threshold = 10**(onset_threshold_dB/20)
    val = torch.max(x_abs)
    val_onset_threshold = val*onset_threshold
    index = (x_abs >= val_onset_threshold).nonzero(as_tuple=False)
    if len(index) == 0:
        return 0
    else:
        first_index = index[0][0]
        return first_index


def get_onset_ch(x,onset_threshold_dB):
    #print(x.shape)
    x_abs = torch.abs(x)
    onset_threshold = 10 ** (onset_threshold_dB / 20)
    val = torch.max(x_abs, dim=1, keepdim=True)[0]  # Max value along each channel
    val_onset_threshold = val * onset_threshold

    #print(val_onset_threshold.shape)
    mask = (x_abs >= val_onset_threshold)
    #print(mask.shape)
    #imgplot = plt.imshow(mask)
    #plt.show()
    
    # Find the indices of non-zero elements along each row (channel)
    non_zero_indices = (mask != 0).nonzero()

    # Initialize a tensor to store the first non-zero indices
    first_non_zero_index = torch.zeros(x.size(0), 1, dtype=torch.float)

    # Loop through each row (channel) to find the first non-zero index
    for i in range(x.size(0)):
        # If there are non-zero elements in the row
        if non_zero_indices.size(0) > 0 and non_zero_indices[:, 0].eq(i).any():
            # Get the first non-zero index along the row
            first_non_zero_index[i] = non_zero_indices[non_zero_indices[:, 0].eq(i)][0][1]
        else:
            # If all elements are zero along the row, set the index to -1
            first_non_zero_index[i] = -1

    return first_non_zero_index.squeeze()

def LPF_sig(x,cutoff_freq,order,fs):
    # Create the LPF filter
    nyquist_freq = fs/2
    b = torch.tensor(signal.firwin(order + 1, cutoff_freq / (nyquist_freq * 2)), dtype=torch.float32)
    b = b.unsqueeze(0).unsqueeze(0)
    
    x_expanded = x.unsqueeze(0).unsqueeze(0)
    filtered_signal_tensor = torch.nn.functional.conv1d(x_expanded, b, bias=None, stride=1, padding=b.shape[-1] // 2).squeeze()
    
    return filtered_signal_tensor

def LPF_sig_ch(x,cutoff_freq,order,fs):
    # Create the LPF filter
    nyquist_freq = fs/2
    b = torch.tensor(signal.firwin(order + 1, cutoff_freq / (nyquist_freq * 2)), dtype=torch.float32)
    b_expanded = b.unsqueeze(0).unsqueeze(0)  # Convert to shape [1, 1, filter_size]

    ## Expand the filter to cover all channels
    #b_expanded = b.expand(1, 1, -1)  # Shape: [num_channels, 1, filter_size]


    x_expanded = x.unsqueeze(1) # add a fake dimantion
    filtered_signal_tensor = torch.nn.functional.conv1d(x_expanded, b_expanded, bias=None, stride=1, padding=b.shape[-1] // 2).squeeze()
    
    return filtered_signal_tensor

def xcorr(x, y):
    """
    Cross-correlation of two 1-D tensors using PyTorch.

    Args:
        x (torch.Tensor): Input sequence 1.
        y (torch.Tensor): Input sequence 2.

    Returns:
        torch.Tensor: Cross-correlation of x and y.
    """

    # Compute cross-correlation using the Fourier domain
    X = torch.fft.fft(x)
    Y = torch.fft.fft(y)
    C = torch.fft.ifft(X * torch.conj(Y))

    # Shift to make cross-correlation result consistent with MATLAB
    C = torch.roll(C, x.shape[0] // 2)

    N = x.shape[0]
    lag_vec = torch.linspace(-N/2, N/2, steps=N)


    return C.real, lag_vec

def get_ITD_CC(x_lr,fs,cutoff_freq,is_plot):
    x_l = x_lr[:,:,0].squeeze()
    x_r = x_lr[:,:,1].squeeze()
    filtered_left  = LPF_sig_ch(x_l,cutoff_freq,10,fs)
    filtered_right = LPF_sig_ch(x_r,cutoff_freq,10,fs)

    

    ITD = torch.zeros(x_lr.shape[0], 1, requires_grad=True)
    #ITD = torch.zeros(x_lr.shape[0], 1)
    for ch_idx in range(filtered_left.shape[0]):
        res, lag_vec = xcorr(filtered_left[ch_idx,:].squeeze(), filtered_right[ch_idx,:].squeeze())
        
        lag_index = torch.argmax(res, dim=0)
        tau_samples = lag_vec[lag_index]
        
        ITD[ch_idx] = 1e6*(tau_samples/fs)
        
        
    if is_plot:
        plt.figure
        plt.plot(ITD)
        plt.xlabel('Incident angle')
        plt.ylabel('Value [micro sec]')
        plt.title('ITD estimation: Cross-correlation')
        plt.grid(True)
        plt.show()
        
    return ITD


def get_ITD_CC_v2(x_lr,fs,cutoff_freq,is_plot):
    x_l = x_lr[:,:,0].squeeze()
    x_r = x_lr[:,:,1].squeeze()
    filtered_left  = LPF_sig_ch(x_l, cutoff_freq, 10, fs)
    filtered_right = LPF_sig_ch(x_r, cutoff_freq, 10, fs)

    itd_list = []
    for ch_idx in range(filtered_left.shape[0]):
        res, lag_vec = xcorr(filtered_left[ch_idx,:].squeeze(), filtered_right[ch_idx,:].squeeze())
        lag_index = torch.argmax(res, dim=0)
        tau_samples = lag_vec[lag_index]
        itd = 1e6 * (tau_samples / fs)
        # Ensure itd is one-dimensional
        itd = itd.view(1) if itd.dim() == 0 else itd
        itd_list.append(itd)

    ITD = torch.cat(itd_list, dim=0).unsqueeze(1).requires_grad_()
    
    if is_plot:
        ITD_tmp = ITD.detach().cpu().numpy()
        title_name = "ITD estimation - Cross correlation method"
        plot_ITD_curve(ITD_tmp,title_name)
        
    return ITD
def get_ITD_onset_threshold(x_lr,fs,trashhold,cutoff_freq,is_plot):
    x_l = x_lr[:,:,0]
    x_r = x_lr[:,:,1]


    filtered_left  = LPF_sig_ch(x_l,cutoff_freq,10,fs)
    filtered_right = LPF_sig_ch(x_r,cutoff_freq,10,fs)

    res_filtered_l = get_onset_ch(filtered_left,trashhold)
    res_filtered_r = get_onset_ch(filtered_right,trashhold)


    ITD = 1e6*((res_filtered_l - res_filtered_r)/fs) # sample diffrence to time [milisec]


   
    if is_plot:
        title_name = "ITD estimation - oneset ditection method"
        plot_ITD_curve(ITD,title_name)
        
        

    return ITD

def get_ITD_group_delay(x_lr,fs,cutoff_freq,is_plot):

    # Calculate the Nyquist frequency
    nyquist_frequency = fs / 2
    
    # Calculate the frequency resolution
    frequency_resolution = fs / x_lr.shape[1]
    
    # Generate the positive frequency vector
    f_vec = torch.arange(0, nyquist_frequency, frequency_resolution)

    
    cutoff_freq_idx = torch.abs(f_vec - cutoff_freq).argmin()

    x_l = x_lr[:,:,0].squeeze().permute(1,0) # [time x ch]
    x_r = x_lr[:,:,1].squeeze().permute(1,0) # [time x ch]

    group_delay_result_L = group_delay_multi_channel(x_l)
    TOA_L = torch.mean(group_delay_result_L[:cutoff_freq_idx,:], 0)
    group_delay_result_R = group_delay_multi_channel(x_r)
    TOA_R = torch.mean(group_delay_result_R[:cutoff_freq_idx,:], 0)

    ITD = (TOA_L - TOA_R) # estimation in micro secs

    if is_plot:
        ITD_tmp = ITD.detach().cpu().numpy()
        ITD_tmp = (ITD_tmp/fs)*1e6
        title_name = "ITD estimation - group delay"
        plot_ITD_curve(ITD_tmp,title_name)
    
    return ITD


def get_muilti_ch_group_delay(x_lr,fs,cutoff_freq,is_plot):
    # Calculate the Nyquist frequency
    nyquist_frequency = fs / 2
    
    # Calculate the frequency resolution
    frequency_resolution = fs / x_lr.shape[1]
    
    # Generate the positive frequency vector
    f_vec = torch.arange(0, nyquist_frequency, frequency_resolution)

    
    cutoff_freq_idx = torch.abs(f_vec - cutoff_freq).argmin()

    x_l = x_lr[:,:,0].squeeze().permute(1,0) # [time x ch]
    x_r = x_lr[:,:,1].squeeze().permute(1,0) # [time x ch]

    group_delay_result_L = group_delay_multi_channel(x_l)[:cutoff_freq_idx,:]
    group_delay_result_R = group_delay_multi_channel(x_r)[:cutoff_freq_idx,:]


    GD = torch.stack((group_delay_result_L,group_delay_result_R),dim=2)
    f_vec_out = f_vec[:cutoff_freq_idx]

    if is_plot:
        GD_tmp = GD.detach().cpu().numpy()
        title_name = "ITD estimation - group delay"
        plt.figure
        plt.semilogx(f_vec_out,GD_tmp[:,:,0])
        plt.grid(True)
        plt.ylim(0,600)
        plt.title(title_name)
        plt.show()
       
    return GD

def plot_ITD_curve(ITD,title_name):
    plt.figure
    plt.plot(ITD)
    plt.xlabel('Incident angle')
    plt.ylabel('Value [micro sec]')
    plt.title(title_name)
    plt.grid(True)
    plt.show()


def group_delay_multi_channel(impulse_response):
    # Get the shape of the input tensor
    N, ch = impulse_response.shape

    # Compute the DFT of the impulse response for all channels
    dft_h = torch.fft.fft(impulse_response, dim=0, norm='forward')

    # Construct the nh[n] sequence for all channels
    n = torch.arange(N).reshape(-1, 1)
    nh = n * impulse_response

    # Compute the DFT of nh[n] for all channels
    dft_nh = torch.fft.fft(nh, dim=0, norm='forward')

    # Calculate the group delay for all channels
    group_delays = torch.real(dft_nh / dft_h)

    # Only retain the positive frequency components
    group_delays = group_delays[:int(N/2)+1, :]

    return group_delays