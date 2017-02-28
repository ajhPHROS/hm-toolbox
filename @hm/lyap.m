function X = lyap(A, C, debug)
%LYAP Solve the Lyapounov equation AX + XA' + C = 0

N = 32;

[x, w] = hm_legpts(N);

if ~exist('debug', 'var')
	debug = false;
end

% Acceleration parameter. 
L = 100;

nrm = norm(A);

% Perform the change of variables and weights
ww = (2 * L * w) .* ( sin(x) ./ (1 - cos(x)).^2 );
xx = L * cot(x ./ 2).^2;

nevals = 0;

X = partialSum(xx, ww, @f);

	function X = partialSum(xx, ww, ef)
		%PARTIALSUM Perform the sum of the function EF at some points. 
		%
		% X = PARTIALSUM(XX, WW, EF) evaluates EF at the points in XX, and
		% combines the evaluations using the weights in WW. 
		if length(xx) == 1
			X = ww * ef(xx);
		else
			mp = ceil(length(xx) / 2);
			
			X = partialSum(xx(1:mp), ww(1:mp), ef) + ...
				partialSum(xx(mp+1:end), ww(mp+1:end), ef);
		end			
	end

	function Y = f(x)	  
	  eA = expm(-x*A, 'pade', 8, nrm * abs(x));
	  nevals = nevals + 1;
	  Y = -eA * C * eA';
	  
	  if debug
		fprintf (' :: %d out of %d evaluations performed\n', nevals, N);
	  end
	end

end

