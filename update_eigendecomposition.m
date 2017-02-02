function [ uu, ss ] = update_eigendecomposition( shift, uu, ss, xx, eta )

% matrix = uu * ( diag( ss ) - shift ) * uu' + shift * eye
% output = input + eta * xx * xx'

dd = size( uu, 1 );
ll = size( uu, 2 );

xx_projected = uu' * xx;    % in the uu basis
xx_leftover = xx - uu * xx_projected;

xx_leftover_norm = norm( xx_leftover );

if ( ( ll < dd ) && ( xx_leftover_norm > 0 ) )

	uu_updated = [ uu, xx_leftover / xx_leftover_norm ];
	ss_updated = [ diag( ss ) + eta * (xx_projected * xx_projected'), eta * xx_projected * xx_leftover_norm ; eta * xx_projected' * xx_leftover_norm, shift + eta * xx_leftover_norm ^ 2 ];

	% svd deals better with repeated eigenvalues than eig
	% however, for negative eigenvalues, the left and right singular vectors may have different signs (so we multiply the singular values by left * right')
	% we may have negative eigenvalues since Warmuth's algorithm works in the log domain
	[ left, values, right ] = svd( ss_updated );
	uu = uu_updated * left;
	ss = diag( values ) .* diag( left * right' );

else

	uu_updated = uu;
	ss_updated = diag( ss ) + eta * (xx_projected * xx_projected');

	% svd deals better with repeated eigenvalues than eig
	% however, for negative eigenvalues, the left and right singular vectors may have different signs (so we multiply the singular values by left * right')
	% we may have negative eigenvalues since Warmuth's algorithm works in the log domain
	[ left, values, right ] = svd( ss_updated );
	uu = uu_updated * left;
	ss = diag( values ) .* diag( left * right' );

end
