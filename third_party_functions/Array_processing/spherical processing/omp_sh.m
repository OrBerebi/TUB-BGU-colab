function [omega, err, x] = omp_sh(anm, opts)
arguments
    anm (:,1) double
    opts.Kmax (1,1) double = inf
    opts.tol (1,1) double = -inf
    opts.omega_grid (:,2) double = sampling_schemes.fliege_maier(29)
    opts.Y_grid (:,:) double = []
    opts.plotFlag (1,1) logical = false
    opts.omega_exp (:,2) double = []
    opts.isComplex (1,1) logical = true
    opts.xRatioThresh (1,1) double = -inf
end
% Author: Tom Shlomo, ACLab BGU, 2020


N = sqrt(size(anm,1))-1;
if isempty(opts.Y_grid)
    opts.Y_grid = shmat(N, opts.omega_grid, opts.isComplex);
end
Q = size(anm,1);
Kmax = min(opts.Kmax, Q);
omega = nan(Kmax,2);
Yh = zeros(Q, Kmax);
r = anm;
err = zeros(Kmax,1);
x = zeros(Kmax);
anm_norm = vecnorm(anm,2,1);
for k=1:Kmax
    [~, omega(k,:)] = sphere_max_abs(r, "omegaGrid", opts.omega_grid, "Ygrid", opts.Y_grid);
    
    if opts.plotFlag
        figure("name", "k="+k);
        h = hammer.surf([],r);
        hold(h.Parent, 'on');
        hammer.plot(omega_exp, [], 'm+');
        hammer.plot(omega(k,:), [], 'rx');
        if k>1
            hammer.plot(omega(1:k-1,:), [], 'r.');
        end
        title("$||r||=$" + norm(r)/norm(anm));
    end
    Yh(:,k) = conj(shmat(N, omega(k,:), true, true));
    
    x(1:k, k) = pinv(Yh(:,1:k))*anm;
    r = anm - Yh(:,1:k)*x(1:k, k);
    err(k) = vecnorm(r,2,1)/anm_norm;
%     fprintf("k=%d: r=%.2f\n", k, err(k));
    if err(k)<=opts.tol
        break
    end
    if k>1 && abs(x(k,k))/abs(x(k-1,k)) < opts.xRatioThresh
        k = k-1; %#ok<FXSET>
        break
    end
end
omega = omega(1:k,:);
err = err(1:k);
x = x(1:k, 1:k);
end

