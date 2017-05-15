function SS = msgproject(S,k)
%msgproject projects to the convex set of constraints sum(S)<=1, S>=0, S<=1/k
%         using binary search.
%
%   S   the vector of eigen-values of the (scaled) projection matrix M
%   k   PCA dimension (desired dimension)

ls=-k; rs=k;
shf=(ls+rs)/2; % set shf = 0.
SS=S+shf; SS(SS<0)=0; SS(SS>1/k)=1/k; % see if capped S has sum leq 1
if sum(SS)<=1 % we are good
    return;
else % sum(SS) > 1
    rs=shf; % search (-k,0) interval
    while abs(ls-rs)>1e-7
        shf=(ls+rs)/2;
        SS=S+shf; SS(SS<0)=0; SS(SS>1/k)=1/k;
        if sum(SS)<=1 % need to search (shf, rs).
            ls=shf;
        else % need to search (ls, shf).
            rs=shf;
        end
    end
end
