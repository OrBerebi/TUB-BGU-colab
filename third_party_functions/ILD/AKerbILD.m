% [err, f_c, n] = AKerbILD(x_l, x_r, f_lim, fs)
%
% Returns the energetic error between a (multi channel) impulse
% response [x] and target function(s) [target] calculated in auditory
% filters of equivalent rectangular bandwidth (ERB).
% x and target must have same length.
%
% For calculating the energetic error, the power spectra of x and the 
% target are filtered with an ERB filterbank. The algorithm integrates over
% each band within the frequency range of the target function and then
% calculates the magnitude difference in dB between x and the target
% function for each integration result.
%
% [Solomons AM (1995) Coloration and binaural decoloration of sound due
% to reflections; Moore BCJ (1995) Hearing]
%
% AKerbError is used by AKpHeadphone. See AKplotDemo.m
%
%
% I N P U T
% x_l         - left ear impulse response 
% x_r       - right ear impulse response 
% f_lim     - upper and lower frequency bounds in Hz.
%             (Default = [50 20000])
% fs        - sampling freuqnecy in Hz. (Default = 44100);
%
% O U T P U T
% err       - error in dB for each ERB filter
% f_c       - mid-frequencies of ERB-Filters
% n         - number of ERB filters used for error analysis
%
% see also: MakeERBFilters, ERBFilterBank (Auditory Toolbox, Slaney 1993)
%
% 08/2008 - initial dev. Zora Schaerer, zora.schaerer@tu-berlin.de
% 07/2013 - minor improvements. Fabian.Brinkmann@tu-berlin.de
% 06/2016 - added frequency dependend filter length for speed up
%           fabian.brinkmann@tu-berlin.de

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function [ILD,f_c,n] = AKerbILD(x_l,x_r, f_lim, fs)

% check if auditory toolbox is in the Matlab search path
%AKdependencies('auditory')

% -------------------------------------------------------- 0. input parsing
% default parameter
if ~exist('f_lim', 'var')
    f_lim = [50 20000];
end
if ~exist('fs', 'var')
    fs = 44100;
end

% throw a warning if input appears to be zero phase
% (all checks if all angles for one spectrum are zero,
%  any checks if this is true for any channel)
xPhase      = any(all(abs(angle(fft(x_l)))      < 1e-10));
targetPhase = any(all(abs(angle(fft(x_r))) < 1e-10));

% don't throw the warning if it is zero phase but a dirac impulse
xDirac      = all( all( x_l(2:end,:) < 1e-10 ) );
targetDirac = all( all( x_r(2:end,:) < 1e-10 ) );

if xPhase && ~xDirac
    warning('AKerbError:Input', 'At least one channel of ''x'' appears to be zero phase. This can produce errors, and can be avoided by circshift(x, [round(size(x,1)/2) 0])')
end
if targetPhase  && ~targetDirac
    warning('AKerbError:Input', 'At least one channel of ''target'' appears to be zero phase. This can produce errors, and can be avoided by circshift(x, [round(size(target,1)/2) 0])')
end


% match IR length
if size(x_l,1)~=size(x_r,1)
    if size(x_l,1)>size(x_r,1)
        x_r(end+1:size(x_l,1),:) = 0;
    else
        x_l(end+1:size(x_r,1),:) = 0;
    end
end

% match IR channels
if size(x_l,2) > size(x_r,2)
    x_r = repmat(x_r, [1 size(x_l,2)]);
end

% get length and channels
N = size(x_l,1);  % ir length
c = size(x_l,2);  % ir num. of channels

% ---------------------------------------- 1. get FIR filter in freq domain
% get ERB filter coefficients
n      = round(21.4*log10(0.004367*f_lim(2)+1)-21.4*log10(0.004367*f_lim(1)+1));    % number of filters
f_c    = ERBSpace(f_lim(1), fs/2, n);                                               % center frequencies
fcoefs = MakeERBFilters(fs,n,f_lim(1));                                             % filter coefficients

% sort with ascending frequency
fcoefs = flipud(fcoefs);
f_c    = flipud(f_c);

% get rid of filters above upper f_lim
n      = find(f_c <= f_lim(2), 1, 'last');
f_c    = f_c(1:n);
fcoefs = fcoefs(1:n,:);

% estimate needed filter length (4 times the largest cycle)
Nmax = 2^nextpow2( 1/f_c(1) * 4 * fs );

% get filter impulse responses
C = zeros(Nmax,1);
C(1) = 1;
C = ERBFilterBank(C',fcoefs)';

% make minimum phase
%C = AKphaseManipulation(C, fs, 'min', 4, 0); %O.B

% filter x with filerbank and calculate error in each band
ILD = zeros(n,c);
C_db_tmp = db(C);
C_db = zeros(size(C_db_tmp));
% ------------------------------------- 2. filter signals and get ERB error
for k = 1:n
    
    % find point where filter decayed for 60 dB
    C_db = C_db_tmp(:,k);
    %C_db = 10*log10(C(:,k)); %O.B
    Nmin = find(flipud(C_db) > max(C_db)-60, 1, 'first');
    Nmin = 2^nextpow2(Nmax - Nmin + 1);
    
    % needed filter length at current center frequency
    Ncur = 2^nextpow2( 1/f_c(k) * 4 * fs );
    Ncur = max([Ncur Nmin N]);
    
    % get indeces of upper and lower frequency limit
    f_lim_id(1) = round(f_lim(1)/(fs/Ncur))+1;
    f_lim_id(1) = max([f_lim_id(1) 2]);         % make sure we dont use the 0 Hz bin
    f_lim_id(2) = round(f_lim(2)/(fs/Ncur))+1;
    
    % filter input signals
    t1 = repmat(abs(fft(C(:,k), Ncur)), [1, c]);
    t2 = (abs(fft(x_l, Ncur)).^2);
    X_L = repmat(abs(fft(C(:,k), Ncur)), [1, c]) .* (abs(fft(x_l, Ncur)).^2);
    X_R = repmat(abs(fft(C(:,k), Ncur)), [1, c]) .* (abs(fft(x_r, Ncur)).^2);
    
    % get energy
    X_L = sum(X_L(f_lim_id(1):f_lim_id(2),:));
    X_R = sum(X_R(f_lim_id(1):f_lim_id(2),:));
    
    % get ERB error
    ILD(k,:) = 10*log10(abs(X_L)./abs(X_R));

end
