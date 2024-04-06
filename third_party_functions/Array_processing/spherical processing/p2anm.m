function anm = p2anm(p, fs, omega_q, r, sphereType, snr_db, N)
arguments
    p (:,:) double
    fs (1,1) double
    omega_q (:,2) double
    r  double
    sphereType (1,1) string {mustBeMember(sphereType, ["rigid", "open", "em32"])}
    snr_db (1,1) double
    N (1,1) double
end
% Author: Tom Shlomo, ACLab BGU, 2020


if sphereType == "em32"
    sphereType = "rigid";
    [th_q, ph_q, ~, r] = sampling_schemes.em32();
    omega_q = [th_q ph_q];
    r = r(1);
end

Q = size(p,2);
% N = floor(sqrt(Q)-1);
c = soundspeed();
pad = round(r/c * fs * 100);
p(end+1:end+1+pad, :) = 0;
NFFT = 2^nextpow2(size(p,1));
P = fft(p, NFFT, 1);
P = P(1:(NFFT/2+1), :);
f = (0 : size(P,1)-1)' * (fs/NFFT);
T = p2anm_mat(f, N, omega_q, r, sphereType, snr_db);
Anm = zeros( size(f,1), (N+1)^2 );
for i=1:size(f,1)
    Anm(i,:) = (T(:,:,i) * P(i,:).').';
end
Anm = sh_freq_complement(Anm, NFFT, 1, 2);
anm = ifft(Anm, [], 1);
anm(size(p,1)+1:end,:) = [];

end

