function [h,bins] = sphere_hist(omega, opts)
arguments
    omega (:,2) double
    opts.tol (1,1) double
    opts.bins (:,2) double
    opts.nbins (1,1) double
    opts.weights (:,1) double = ones(size(omega,1),1)
    opts.type (1,1) string {mustBeMember(opts.type, ["nearest", "overlap"])}
    opts.plotFlag (1,1) logical = false
end
% Author: Tom Shlomo, ACLab BGU, 2020


%% parse inputs
if ~isfield(opts, "type")
    if isfield(opts, "tol")
        opts.type = "overlap";
    else
        opts.type = "nearst";
    end
end

if ~isfield(opts, "bins") || isempty(opts.bins)
    bins = sampling_schemes.golden_spiral(opts.nbins);
else
    bins = opts.bins;
end

%%

m = size(bins,1);
h = zeros(m,1);
switch opts.type
    case "overlap"
        for i=1:m
            d = angle_between(bins(i,:), omega);
            I = d <= opts.tol;
            h(i) = sum(opts.weights(I));
        end
    case "nearest"
        error("not yet implemented");
end 

if opts.plotFlag
    hammer.scatter3(bins, [], h, 1000, h, '.');
    colormap(flipud(gray));
    colorbar;
end

end

