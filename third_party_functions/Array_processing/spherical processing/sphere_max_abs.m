function [g, omega, traj] = sphere_max_abs(fnm, opts)
arguments
    fnm double
    opts.newtonFlag (1,1) logical = true
    opts.newtonTol (1,1) double {mustBePositive, mustBeFinite} = 2*pi*1e-4
    opts.newtonMaxIter (1,1) double {mustBeInteger, mustBePositive} = 10
    opts.isComplexSH (1,1) logical = true
    opts.plotFlag (1,1) logical = false
    opts.parent (1,1) matlab.graphics.axis.Axes
    opts.omegaGrid (:,2) double = sampling_schemes.fliege_maier(29);
    opts.Ygrid double = [];
    opts.normalization (1,1) string {mustBeMember(opts.normalization, ["none", "rho", "directivity", "normalizeddirectivity", "music"])} = "none"
    opts.wbar (1,1) logical = size(fnm,2) > 100
    opts.buffer_size = 1e7;
end
% Author: Tom Shlomo, ACLab BGU, 2020


N = ceil(sqrt(size(fnm,1))-1);
if isempty(opts.Ygrid)
    opts.Ygrid = shmat(N, opts.omegaGrid, opts.isComplexSH, false);
end
buffer_size = round(opts.buffer_size/size(opts.Ygrid, 1));
I0 = 0:(buffer_size-1);
num_columns = size(fnm, 2);
g = zeros(1, num_columns);
f = zeros(1, num_columns);
j = zeros(1, num_columns);
for i=1:buffer_size:num_columns
    I = I0+i;
    if I(end)>num_columns
        I(I>num_columns)=[];
    end
    f_all_dirs = opts.Ygrid*fnm(:,I);
    [g(I),j(I)] = max(abssq(f_all_dirs), [], 1);
    f(I) = f_all_dirs(sub2ind(size(f_all_dirs), j(I), find(I)));
end
omega = opts.omegaGrid(j,:);

if opts.newtonFlag
    if opts.wbar
        H = wbar();
    end
    traj = repmat(struct('omega', [], 'g', []), size(fnm,2), 1);
    for i=1:size(fnm,2)
        [g_traj, omega_traj] = sphere_abs_newton(fnm(:,i), omega(i,:),...
            "f0", f(i), ...
            "isComplex", opts.isComplexSH,...
            "maxIter", opts.newtonMaxIter,...
            "tol", opts.newtonTol,...
            "Y0", opts.Ygrid(j(i),:), ...
            "plotFlag", false, ...
            "stopIfNotNegativeDefinite", true);
        [g(i), k] = max(g_traj);
        omega(i,:) = omega_traj(k,:);
        if opts.wbar && (mod(i,1000)==0 || i==size(fnm,2)); wbar(i, size(fnm,2),H); end
        if nargout==3
            traj(i).omega = omega_traj;
            traj(i).g = g_traj;
        end
    end
end

switch opts.normalization
    case "none"

    case "rho"
        g = sqrt(g*4*pi*(N+1)^(-2))./vecnorm(fnm,2,1);
    case "directivity"
        g = g./( vecnorm(fnm,2,1).^2 );
        warning("hasn't been tested yet");
    case "normalizeddirectivity"
        g = g./( vecnorm(fnm,2,1).^2 *(N+1)^2 );
        warning("hasn't been tested yet");
    case "music"
        g = 1./( 1 - g./( vecnorm(fnm,2,1).^2 *(N+1)^2 ) );
        warning("hasn't been tested yet");
    otherwise
        error("unknown normalization option: %s", opts.normalization);
end

if opts.plotFlag
    if ~isfield(opts, "parent")
        opts.parent = gca();
    end
    hammer.surf([],fnm(:,1),[], @(f) 10*log10(abssq(f)), true, 'Parent', opts.parent);
    hold on;
    hammer.plot(omega(1,1), omega(1,2), 'r.-', 'Parent', opts.parent);
    warning("plot does not reflect normalization");
end

end

