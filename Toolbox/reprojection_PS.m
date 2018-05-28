function SE =  reprojection_PS(S,rho,N,I,imask)

	[nrows,ncols,nimgs] = size(I);
	S1 = reshape(S(:,1),[1 1 nimgs]);
	S2 = reshape(S(:,2),[1 1 nimgs]);
	S3 = reshape(S(:,3),[1 1 nimgs]);

	Ireproj = bsxfun(@times,S1,rho.*N(:,:,1))+bsxfun(@times,S2,rho.*N(:,:,2))+bsxfun(@times,S3,rho.*N(:,:,3));

	I_vect = reshape(I,nrows*ncols,nimgs);
	I_vect = I_vect(imask,:);
	Ireproj_vect = reshape(Ireproj,nrows*ncols,nimgs);
	Ireproj_vect = Ireproj_vect(imask,:);

	SE = zeros(nrows,ncols);
	SE(imask) = 0.5*mean((I_vect-Ireproj_vect).^2,2);
end
