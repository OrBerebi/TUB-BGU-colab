function [th, ph, weights] = t_design(N)
%GETTDESIGN Returns the spherical coordinates of minimal T-designs
%
%   GETTDESIGN returns the unit vectors and the spherical coordinates
%   of t-designs, which constitute uniform arrangements on the sphere for
%   which spherical polynomials up to degree t can be integrated exactly by
%   summation of their values at the points defined by the t-design.
%   Designs for order up to t=21 are stored and returned. Note that for the
%   spherical harmonic transform (SHT) of a function of order N, a spherical
%   t-design of t>=2N should be used (or equivalently N=floor(t/2) ), since
%   the integral evaluates the product of the spherical function with
%   spherical harmonics of up to order N. The spherical coordinates are
%   given in the [azi1 elev1; azi2 elev2; ...; aziQ elevQ] convention.
%
%   The designs have been copied from:
%       http://neilsloane.com/sphdesigns/
%   and should be referenced as:
%       "McLaren's Improved Snub Cube and Other New Spherical Designs in
%       Three Dimensions", R. H. Hardin and N. J. A. Sloane, Discrete and
%       Computational Geometry, 15 (1996), pp. 429-441.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, archontis.politis@aalto.fi, 10/11/2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modified by Tom Shlomo, Dec 2019

t = 2*N;

if t>21
    error('Designs of order greater than 21 are not implemented.')
elseif t<1
    error('Order should be at least 1.')
end

persistent t_designs
if isempty(t_designs)
    load('+sampling_schemes/t_designs_1_21.mat', 't_designs');
end
vecs = t_designs{t};
[th, ph] = c2s(vecs(:,1), vecs(:,2), vecs(:,3));
weights = repmat( 4*pi/length(th) , size(th) );

if nargout==1
    th = [th ph];
end
end