
%% Extracts the best maximal kk-dimensional subspace from a solution
%
% kk - the dimension of the subspace which we seek
% left, right, values - "nontrivial" signular vectors and singular values
%                   of the iterate
%%

function [ U, V ] = top_subspace( k, U, S, V )
[ ~, indices ] = sort( S, 'descend' );
U = U( :, indices( 1:k ) );
V = V( :, indices( 1:k ) );

U  = U  / sqrtm( U'  * U  );
V = V / sqrtm( V' * V );

end
