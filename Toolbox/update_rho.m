function rho = update_rho(I,zu,zv,S,mask,S1,S2,S3)
% Least squares update of the albedo
% In :
%  I : nrows x ncols x nimgs
%  S : nimgs x 3
%  zu and zv : nrows x nimgs
%  mask : nrows x ncols
%  S1,S2,S3 : 1 x 1 x nimgs
% Out :
%  rho : nrows x ncols

	if(nargin < 4)
		disp('not enough input arguments');
		return;
	end

	[nrows,ncols,nimgs] = size(I);
	n = nrows*ncols;

	if(nargin < 5)
		mask = ones(nrows,ncols);
	end
	if(nargin <7)
		S1 = reshape(S(:,1),[1 1 nimgs]); % S(:,1) for bsxfun
		S2 = reshape(S(:,2),[1 1 nimgs]); % S(:,2) for bsxfun
		S3 = reshape(S(:,3),[1 1 nimgs]); % S(:,3) for bsxfun
	end

	% Some quantities
	minus_I_times_S1_dot_zu = I.*bsxfun(@times,S1,-zu); % nrows x ncols x nimgs
	minus_I_times_S2_dot_zv = I.*bsxfun(@times,S2,-zv); % nrows x ncols x nimgs
	I_times_S3 = bsxfun(@times,S3,I); % nrows x ncols x nimgs
	minus_S1_dot_zu = bsxfun(@times,S1,-zu); % nrows x ncols x nimgs
	minus_S2_dot_zv = bsxfun(@times,S2,-zv); % nrows x ncols x nimgs

	% Numerator
	num = sum(minus_I_times_S1_dot_zu+minus_I_times_S2_dot_zv+I_times_S3,3); % nrows x ncols
	denom = sum((bsxfun(@plus,S3,minus_S1_dot_zu+minus_S2_dot_zv)).^2,3);

	% Result
	rho = sqrt(1+zu.^2+zv.^2).*num./denom;
	rho(denom == 0) = 0;
	rho(mask == 0) = 0;
end
