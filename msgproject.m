function SS = msgproject(S,k)
%msgproject projects to the convex set of constraints sum(S)=1, S>=0, S<=1/k
%           using binary search. Note that this algorithm is not efficient.
%           For an efficient algorithm please refer to Arora et al 2013
%
%   S   the vector of eigen-values of the (scaled) projection matrix M
%   d   dimension of the ambient space
%   k   PCA dimension (desired dimension)
%   eps the maximum error that can be tolerated in projection

ls=-k; rs=k;
while abs(ls-rs)>1e-6
    shf=(ls+rs)/2;
    SS=S+shf; SS(SS<0)=0; SS(SS>1/k)=1/k;
    if sum(SS)>1
        rs=shf;
    else
        ls=shf;
    end
end
end
