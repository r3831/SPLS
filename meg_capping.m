function [ uu, ss ] = meg_capping( uu, ss, dd, kk )

bound = 1 / ( dd - kk );

while true

	ll = length( ss );

	% normalize the eigenvalues
	ss = ss * ( 1 - ( dd - ll ) * bound ) / sum( ss );

	% check if any exceed the bound
	[ value, index ] = max( ss );
	if ( value < bound ); break; end

	% set the maximum eigenvalue to bound, and loop
	% recall that any directions which are not in the basis are assumed to have an eigenvalue of bound
	if ( index == 1 )
		uu = uu( :, 2:ll );
		ss = ss( 2:ll  );
	elseif ( index == ll )
		uu = uu( :, 1:( ll - 1 ) );
		ss = ss( 1:( ll - 1 ) );
	else
		uu = [ uu( :, 1:( index - 1 ) ), uu( :, ( index + 1 ):ll ) ];
		ss = [ ss( 1:( index - 1 )  ) ; ss( ( index + 1 ):ll ) ];
	end

end

ss( ss < eps ) = eps;
