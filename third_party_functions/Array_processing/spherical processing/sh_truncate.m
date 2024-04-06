function fnm = sh_truncate(fnm, N, dim, padFlag)

% Author: Tom Shlomo, ACLab BGU, 2020

if nargin<4
    padFlag = 0;
end
if nargin<3 || isempty(dim)
    dim = find(size(fnm)>1,1);
    if isempty(dim)
        dim = 1;
    end
end
Np1sq = (N+1)^2;
if size(fnm,dim)==Np1sq
    return
elseif size(fnm,dim)<Np1sq
    if padFlag
        i = size(fnm,dim)+1 : Np1sq;
        S.subs = repmat({':'},1,ndims(fnm));
        S.subs{dim} = i; % the third row
        S.type = '()';
        fnm = subsasgn(fnm, S, 0 );
        return
    else
        error('impossible');
    end
end
i = 1:Np1sq;
S.subs = repmat({':'},1,ndims(fnm));
S.subs{dim} = i; % the third row
S.type = '()';
fnm = subsref(fnm,S);

end

