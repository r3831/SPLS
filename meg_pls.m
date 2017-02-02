function [ U, S, V ] = meg_pls( k, U, S, V, X, Y, eta )

W = [ U ; V ];

dx = size( U,  1 );
d = size( W, 1 );

bound = 1 / ( d - k );
logbound = log( bound );

% update matrix = exp( log( matrix ) - eta * [ 0, Y * X' ; X * Y', 0 ] )
S = log( S );
[ W, S ] = update_eigendecomposition( logbound, W, S, [ X ;  Y ], -0.5 * eta );
[ W, S ] = update_eigendecomposition( logbound, W, S, [ X ; -Y ],  0.5 * eta );
S = exp( S );

% cap the maximum eigenvalue to bound
[ W, S ] = meg_capping( W, S, d, k );

U  = W(            1:dx, : );
V = W( ( dx + 1 ):d,   : );
