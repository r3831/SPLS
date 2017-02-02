function [ U, V ] = meg_pls_solution( keff, U, dS, V )

[ ~, indices ] = sort( dS );

U  = U(  :, indices( 1:keff ) );
V = V( :, indices( 1:keff ) );

U  = U  / sqrtm( U'  * U  );
V = V / sqrtm( V' * V );
