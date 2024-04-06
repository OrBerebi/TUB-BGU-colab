function p2anm = p2anm_mat(f, N, omega_q, r, sphereType, snr_db, Sigma_n, Sigma_anm)
arguments
    f (:,1) double % Hz
    N (1,1) double
    omega_q (:,2) double
    r double
    sphereType (1,1) string {mustBeMember(sphereType, ["rigid", "open", "em32"])}
    snr_db (1,1) double = []
    Sigma_n (:,:) double = 0;
    Sigma_anm (:,:) double = inf;
end
% Author: Tom Shlomo, ACLab BGU, 2020



if sphereType == "em32"
    sphereType = "rigid";
    [th_q, ph_q, ~, r] = sampling_schemes.em32();
    omega_q = [th_q ph_q];
    r = r(1);
end
Q = size(omega_q,1);
F = size(f,1);
kr = 2*pi*f*r/soundspeed();
if ~isempty(snr_db)
    Sigma_anm = 1;
    snr = 10^(snr_db/10);
    Sigma_n = Sigma_anm*4*pi/snr;
end
if isscalar(Sigma_anm)
    Nmax = max([ceil(max(kr))+1, N]);
    Sigma_anm = eye((Nmax+1)^2)*Sigma_anm;
else
    Nmax = size(Sigma_anm,1);
end
if isscalar(Sigma_n)
    Sigma_n = Sigma_n*eye(Q);
end

% if isempty(Sigma_anm) && (N+1)^2 > Q
%     warning("Under determined system, consider using a bayesian estimator instead");
% end
% if isempty(Sigma_anm) && max(kr) > N
%     warning("U
% end
% if isscalar(Sigma_n)
%     Sigma_n = Sigma_n*eye(Q);
% end
% if isscalar(Sigma_anm)
%     Sigma_anm = Sigma_amm*eye((N+1)^2);
% end
Y = shmat(Nmax, omega_q);
B = bn(Nmax, kr, "sphereType", sphereType);
% anm2p = Y .* permute(B, [3 2 1]);
p2anm = zeros( (N+1)^2, Q, F );
if isfinite(trace(Sigma_anm))
    Sigma_anm_trunc = sh_truncate(Sigma_anm, N, 1);
    for i=1:F
        H = Y.* B(i,:);
        
        p2anm(:,:,i) = ( Sigma_anm_trunc * H' )/ (H*Sigma_anm*H' + Sigma_n);
    end
elseif trace(Sigma_n)>0
    Sigma_n_inv = inv(Sigma_n);
    Sigma_anm_inv = inv(Sigma_anm);
    
    for i=1:F
        H = Y.* B(i,:);
        Hh_trunc = sh_truncate(H', N, 1);
        p2anm(:,:,i) = (Hh_trunc*Sigma_n_inv*H + Sigma_anm_inv)\H'*Sigma_n_inv; %#ok<MINV>
    end
else
    if ~((N+1)^2 >= Q && max(kr) <= N )
        warning("Spatial aliasing/under determined system");
    end
    for i=1:F
        H = Y.* B(i,:);
        p2anm(:,:,i) = sh_truncate(pinv( H ), N, 1);
    end
end

% p2anm = p2anm( 1:(N+1)^2, :, :)';



end

