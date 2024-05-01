% Discription: Load the HRTF in sofa format, convert to Earo format,
% calculate N'th order SH coefficients

% Written by: O.B 04 Apr 2024

clc
clear
close all
restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED
set(0,'DefaultFigureWindowStyle','docked')
addpath(genpath("/Users/orberebi/Documents/GitHub/TUB-BGU-colab"));

% set graphics defaults:
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultLineMarkerSize', 10);
set(0, 'defaultAxesXGrid', 'on')
set(0, 'defaultAxesYGrid', 'on')
set(0, 'defaultAxesZGrid', 'on')
set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(groot,'defaultAxesZMinorGrid','on','defaultAxesZMinorGridMode','manual');
set(0,'DefaultTextFontSize',20) % Sets font size
set(0,'defaultAxesFontSize',20)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1);


% Flags
is_plot = true;
is_save = false;


% HRTF file path and name
file_path_HRTF = "/Users/orberebi/Documents/GitHub/TUB-BGU-colab/HRTF/KU100_HRIR_L2702.sofa";

% N parameters
N_low  = 1;
N_high = 42;

% Cutoff frequency for MaleS preprocessing
cutOffFreq = 1.5e3;
% Frequency band for evaluating ILD curves
f_band = [1.5e3, 20e3];

% Omega_az
ph_DOA_omega_az = deg2rad((-180:180)).';
th_DOA_omega_az = deg2rad(90 * ones(size(ph_DOA_omega_az)));
omega_az  = [th_DOA_omega_az, ph_DOA_omega_az];
clear ph_DOA_res ph_DOA_max ph_DOA_omega_az th_DOA_omega_az


% Save data path
save_name = "/Users/orberebi/Documents/GitHub/TUB-BGU-colab/matlab_saved_data/25_04_24/";
mkdir(save_name)
[fPath, fName, fExt] = fileparts(file_path_HRTF);
tmp = "/"+fName + "_N"+num2str(N_low)+".mat";
save_name = save_name+tmp;
clear tmp


disp("Convert .sofa to earo object")
hobj = sofaToEaro(convertStringsToChars(file_path_HRTF));
fs = hobj.fs;
nfft = hobj.taps;
nfft = nfft*4;
hobj.taps = nfft;
f_vec = linspace(0,fs/2,nfft/2+1).';
f_vec(1) = f_vec(2)/2;

% Get the angles for the sphirical grid
omega = [hobj.sourceGrid.elevation,hobj.sourceGrid.azimuth];




disp("Calculating high order SH coefficients")
Hnm_high = get_SH_coeff(hobj,N_high,nfft); % SH x freq x ears

disp("Calculating high order Y matrices")
Y_high_lebedev = sh2(N_high,omega(:,1),omega(:,2)).';   % SH x Omega
Y_high_az      = sh2(N_high,omega_az(:,1),omega_az(:,2)).';                         % SH x Omega_az

disp("Calculating high order ear pressure signals")
[p_f_high_lebedev,p_t_high_lebedev] = get_p_signals(Hnm_high,Y_high_lebedev,nfft);
[p_f_high_az,p_t_high_az]           = get_p_signals(Hnm_high,Y_high_az,nfft);

disp("Calculating low order SH coefficients")
%Hnm_low = get_SH_coeff(hobj,N_low,nfft,Y_high_lebedev,cutOffFreq,f_vec); % SH x freq x ears
Hnm_low = get_SH_coeff(hobj,N_low,nfft); % SH x freq x ears

disp("Calculating low order Y matrices")
Y_low_lebedev = sh2(N_low,omega(:,1),omega(:,2)).';   % SH x Omega
Y_low_az      = sh2(N_low,omega_az(:,1),omega_az(:,2)).';                         % SH x Omega_az

disp("Calculating low order ear pressure signals")
[p_f_low_lebedev,p_t_low_lebedev] = get_p_signals(Hnm_low,Y_low_lebedev,nfft);
[p_f_low_az,p_t_low_az]           = get_p_signals(Hnm_low,Y_low_az,nfft);


disp("Calculating low order SH MagLS coefficients")
Hnm_mls = get_MagLS_coeff(hobj,nfft, N_low, Y_low_lebedev, Y_high_lebedev, cutOffFreq);
disp("Calculating low order MagLS ear pressure signals")
[p_f_mls_lebedev,p_t_mls_lebedev] = get_p_signals(Hnm_mls,Y_low_lebedev,nfft,true);
[p_f_mls_az,p_t_mls_az]           = get_p_signals(Hnm_mls,Y_low_az,nfft,false);


disp("calculate the NMSE error")
e_mse_mls = mean(squeeze(mean(10*log10((abs(p_f_high_lebedev-p_f_mls_lebedev).^2)./abs(p_f_high_lebedev).^2),1)),2);
e_mse_ls  = mean(squeeze(mean(10*log10((abs(p_f_high_lebedev-p_f_low_lebedev).^2)./abs(p_f_high_lebedev).^2),1)),2);

disp("calculate the Magnitude error")
e_mag_mls = mean(squeeze(mean(10*log10((abs(abs(p_f_high_lebedev)-abs(p_f_mls_lebedev)).^2)./abs(p_f_high_lebedev).^2),1)),2);
e_mag_ls  = mean(squeeze(mean(10*log10((abs(abs(p_f_high_lebedev)-abs(p_f_low_lebedev)).^2)./abs(p_f_high_lebedev).^2),1)),2);

disp("calculate the Magnitude error")
[ILD_ref,f_c] = get_ILD_curves(p_t_high_az,f_band,fs);
[ILD_mls,~] = get_ILD_curves(p_t_mls_az,f_band,fs);
[ILD_ls,~]  = get_ILD_curves(p_t_low_az,f_band,fs);
%[X_ref,f_c] = get_gammatone(p_t_high_az,f_band,fs);

e_ILD_ls = squeeze(mean(abs(ILD_ref - ILD_ls),1)).';
e_ILD_mls = squeeze(mean(abs(ILD_ref - ILD_mls),1)).';
e_ILD_ls_freq  = squeeze(mean(abs(ILD_ref - ILD_ls),2));
e_ILD_mls_freq = squeeze(mean(abs(ILD_ref - ILD_mls),2));


ang_vec = rad2deg(omega_az(:,2));

if is_save
    save(save_name,"p_f_high_lebedev","ILD_ref","f_band","Hnm_low","Hnm_mls","Hnm_high","fs","N_low","N_high","f_vec","ang_vec","Y_high_lebedev","Y_high_az","Y_low_lebedev","Y_low_az","omega","omega_az");

end


if is_plot
    disp("Plotting the HRIR for a certain direction")
    t = linspace(0,nfft/fs,nfft);
    

    figure
    f_vec = linspace(0,fs/2,nfft/2+1).';
    f_vec(1) = f_vec(2)/2;
    semilogx(f_vec,e_mse_ls(:,1));
    hold on
    semilogx(f_vec,e_mse_mls(:,1));
    semilogx(f_vec,e_mag_ls(:,1));
    semilogx(f_vec,e_mag_mls(:,1));
    legend(["NMSE - LS","NMSE - MagLS","Magnitude Error - LS","Magnitude Error - MagLS"],"Location","northwest")
    grid on
    xlim([50 20e3])
    title("Frequency analysis averaged over dense lebedev grid and over both ears")
    xlabel("Frequency [Hz]")
    ylabel("Error [dB]")


    figure
    
    plot(ang_vec,e_ILD_ls)
    hold on
    plot(ang_vec,e_ILD_mls)
    legend(["LS","MagLS"])
    grid on
    xlim([-180 180])
    ylim([0 12])
    title("ILD error")
    xlabel("Incident Angle [deg]")
    ylabel("Error [dB]")



    figure
    semilogx(f_c,e_ILD_ls_freq,'--o');
    hold on
    semilogx(f_c,e_ILD_mls_freq,'--o');
    legend(["LS","MagLS"],"Location","northwest")
    grid on
    xlim(f_band)
    ylim([0 12])
    title("ILD error frequency analysis averaged over Omega az grid")
    xlabel("Frequency [Hz]")
    ylabel("Error [dB]")

    
end


function Hnm = get_SH_coeff(hobj,N,nfft,Y_high,cutOffFreq,f_vec)
    if nargin<4
        Y_high = [];
        cutOffFreq = 0;
    end

    hobj.shutUp = true;
    if strcmp(hobj.dataDomain{1}, 'TIME')
        hobj = hobj.toFreq(nfft);
        hobj.data = hobj.data(:,1:end/2+1,:);
        hobj.taps = nfft;
    end
    
    if strcmp(hobj.dataDomain{2}, 'SPACE')
        H = double(hobj.data);
        degrees = size(H,1);
        [gridData_Inf, ~, ~]    = sofia_lebedev(degrees,0);

        ph = gridData_Inf(:,1); 
        th = gridData_Inf(:,2);
        a  = gridData_Inf(:,3)*4*pi;

        Y = shmat(N,[th,ph]);   % SH x Omega
        Yp = Y'*diag(a);

        
        for nn=1:size(H,3)       % iterate through mics
            %tmp(:,:,nn) = (double(squeeze(H(:,:,nn)).')*Yp).';
            tmp(:,:,nn) = Yp*H(:,:,nn);
        end
        Hnm = tmp;


        hobj = hobj.toSH(N,'SRC');

    end

    Hnm = hobj.data; 
    if N < 3 && ~isempty(Y_high)
        k = 1;
        fc_high = cutOffFreq*2^(k/2);
        
        % find cut-off freq index
        srcharray = abs(f_vec-fc_high);
        cutOffFreqInd_high = find(srcharray==min(srcharray));
        [tmp_l,tmp_r] = diff_field_eq(Y_high.',H,Hnm(:,:,1).',Hnm(:,:,2).',cutOffFreqInd_high);
        Hnm(:,:,1) = tmp_l.';
        Hnm(:,:,2) = tmp_r.';
    end
end
function Hnm = get_MagLS_coeff(hobj,nfft, N, Y_N, Y_high, cutOffFreq)
    hobj.shutUp = true;
    if strcmp(hobj.dataDomain{1}, 'TIME')
        hobj = hobj.toFreq(nfft);
        hobj.data = hobj.data(:,1:end/2+1,:);
        hobj.taps = nfft;
    end
    H = hobj.data;
    fs = hobj.fs;
    f_vec = linspace(0,fs/2,size(H,2)).';
    Hnm = [];
    [H_l_nm_MagLS, H_r_nm_MagLS] = computeMagLS_imp(H, f_vec, N, Y_N.', cutOffFreq, Y_high.');
    Hnm = cat(3,Hnm,H_l_nm_MagLS);
    Hnm = cat(3,Hnm,H_r_nm_MagLS);
    Hnm = permute(Hnm,[2,1,3]);
end
function [p_f,p_t] = get_p_signals(hnm,Y,nfft,flag_circ)
    if nargin<4
        flag_circ = false;
    end
    p_f = [];
    p_hat_omega_l = Y*hnm(:,:,1);  
    p_hat_omega_r = Y*hnm(:,:,2);
    p_f = cat(3,p_f,p_hat_omega_l);
    p_f = cat(3,p_f,p_hat_omega_r);
    p_f_tmp = p_f;
    p_f_tmp(:, end+1 : nfft, :) = 0;
    p_t = ifft(p_f_tmp, [], 2, 'symmetric');
    
    if flag_circ
        p_t  = circshift(p_t, floor(nfft / 4), 2);
    end

end
function [ILD_no_av,f_c] = get_ILD_curves(P_t,f_band,fs)
    Pl = P_t(:,:,1).';
    Pr = P_t(:,:,2).';
    [ILD_no_av,f_c,~] = AKerbILD(Pl, Pr, f_band, fs);
end

function [x_pos_grad,f_c] = get_gammatone(P_t,f_band,fs)
    Pl = P_t(:,:,1).';
    Pr = P_t(:,:,2).';
    [x_pos_grad,f_c] = Gammatone_filt(Pl, Pr, f_band, fs);
end