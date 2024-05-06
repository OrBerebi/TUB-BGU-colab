function [H_l_nm_MagLS, H_r_nm_MagLS] = computeMagLS_imp_phase(H, f_vec, N, Y_N, cutOffFreq, Y_high)
% function computeMagLS computes the MagLS representation of the HRTFs
%
% Implementation based on Ambisonics book (https://link.springer.com/book/10.1007/978-3-030-17207-7)
%   Eqs.(4.57)-(4.59), for N<3 until (4.63)
%
% Inputs:
%       H          : HRTF in space/frequncy domain [# of directions X # of frequency bins X 2 (left/right)]
%       f_vec      : Frequencies vector [Hz]
%       N          : SH order
%       Y_N        : Complex normalized SH functions matrix of order N [(N+1)^2 X # of directions ]
%       cutOffFreq : Frequency cut-off for MagLS [Hz] (set at 2KHz by deafult)
%       Y_high     : Complex normalized SH functions matrix of high order M>>N (used for Diffuse-Field Covariance Constraint, when N<3)
% Output:
%       H_l_nm_MagLS : SH coefficient of left HRTF [# of frequency bins X (N+1)^2]
%       H_r_nm_MagLS : SH coefficient of right HRTF [# of frequency bins X (N+1)^2]
%
% Written by Zamir Ben-Hur
% Nov 2020
%
% Modified by Or B May 2024
% Added Christoph Hold phase Contiunuation from 2023 DAGA

% initials
if nargin < 6 % don't use Diffuse-Field Covariance Constraint
    Y_high = [];
    flag = false;
end
% H = H(get_90,:,:);
% Y_high = Y_high(:,get_90);
% Y_N = Y_N(:,get_90);
fs = f_vec(end)*2;
hf_delay = [0,0];

k = 1;
fc_low = cutOffFreq*2^(-k/2);
fc_high = cutOffFreq*2^(k/2);

% find cut-off freq index
srcharray = abs(f_vec-fc_low);
cutOffFreqInd_low = find(srcharray==min(srcharray));
srcharray = abs(f_vec-fc_high);
cutOffFreqInd_high = find(srcharray==min(srcharray));

srcharray = abs(f_vec-cutOffFreq);
cutOffFreqInd = find(srcharray==min(srcharray));


pY = pinv(Y_N); 

H_l = squeeze(double(H(:,:,1))).'; % left HRTF
H_r = squeeze(double(H(:,:,2))).'; % right HRTF

H_l_nm_MagLS = zeros(size(H,2),(N+1)^2);
H_r_nm_MagLS = zeros(size(H,2),(N+1)^2);

phi_l_mod = zeros(size(H_l));
phi_r_mod = zeros(size(H_r));


% Compute MagLS

% for freq < cutOffFreq, standard SH transform (Eq. 4.57)
H_l_nm_MagLS(1:cutOffFreqInd,:) = H_l(1:cutOffFreqInd,:)*pY; % SH transform - left ear
H_r_nm_MagLS(1:cutOffFreqInd,:) = H_r(1:cutOffFreqInd,:)*pY; % SH transform - right ear


phi_l_mod(1:cutOffFreqInd,:) = angle(H_l_nm_MagLS(1:cutOffFreqInd,:) * Y_N );
phi_r_mod(1:cutOffFreqInd,:) = angle(H_r_nm_MagLS(1:cutOffFreqInd,:) * Y_N );

% get the delta (from prediction frame)
% for each angle
delta_phi_l = mean(diff(unwrap(phi_l_mod(1:cutOffFreqInd,:)), 1), 1);
delta_phi_r = mean(diff(unwrap(phi_r_mod(1:cutOffFreqInd,:)), 1), 1);

% apply (additional) group hf delay
delta_w = 2*pi*f_vec(2);
delta_phi_l = (delta_phi_l - delta_w * (hf_delay(1) / fs)).';
delta_phi_r = (delta_phi_r - delta_w * (hf_delay(2) / fs)).';


% for freq >= cutOffFreq, MagLS (Eqs. 4.58, 4.59)

for freqInd = cutOffFreqInd + 1 : size(H,2)
    
    % left ear
    phi_est_l = angle(Y_N.'*H_l_nm_MagLS(freqInd-1,:).'); % Eq. (4.58)
    phi_est_l = phi_est_l + delta_phi_l;
    H_l_nm_MagLS(freqInd,:) = (abs(H_l(freqInd,:)) .* exp(1i*phi_est_l.')) * pY; %  Eq. (4.59)
    
    % right ear
    phi_est_r = angle(Y_N.'*H_r_nm_MagLS(freqInd-1,:).'); % Eq. (4.58)
    phi_est_r = phi_est_r + delta_phi_r;
    H_r_nm_MagLS(freqInd,:) = (abs(H_r(freqInd,:)) .* exp(1i*phi_est_r.')) * pY; %  Eq. (4.59)
end

%for N<3 use Diffuse-Field Covariance Constraint
flag =false;
if flag
    if N < 3 && ~isempty(Y_high)
        [H_l_nm_MagLS,H_r_nm_MagLS] = diff_field_eq(Y_high,H,H_l_nm_MagLS,H_r_nm_MagLS,cutOffFreqInd_high);


end
end
end