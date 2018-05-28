function [func_f,zu,zv,grad_f] = func_f(z,I,rho,M,imask,S,II1,II2,JJ1,JJ2,Ivect,grad_approx)
    
	nimgs = size(I,3);
	[func_f,A,b,dz_curr,zu,zv] = func_Ab(z(imask),I,rho,M,imask,S,II1,II2,JJ1,JJ2,Ivect);

    % Compute gradient
    if(nargout>=4)
        npix = length(imask);
        
        D = A*M;
        
        if grad_approx
            grad_f = D'*(A*M*z(imask)-b)./(nimgs);
        else
            D_t = D';
            rho_vec = rho(imask);
            M_t = M';
            zgrad = [zu(imask), zv(imask)];
            zgrad_t = zgrad';
            MMz = bsxfun(@times,M_t,(zgrad_t(:))');
            
            MMzSum = MMz(:,1:2:end-1)+MMz(:,2:2:end);
            
            grad = [zgrad, -ones(npix,1)];
            
            tmp1 = (ones(nimgs,1)*(rho_vec./(dz_curr(imask)).^3)').*(S*grad');
            
            tmp2 = cell(npix,1);
            tmp3 = cell(npix,1);
            tmp4 = cell(npix,1);
            for j=1:npix
                tmp2{j,1} = tmp1(:,j)';
            end
            for j=1:npix
                tmp3{j,1} = MMzSum(:,j);
            end
            parfor j=1:npix
                tmp4{j,1} = sparse(tmp3{j,1}*tmp2{j,1});
            end
            D_t = D_t+cell2mat(tmp4');
            
            grad_f = D_t*(A*M*z(imask)-b)./(nimgs);
        end
    end
end

function [func_f,A,b,dz,zu,zv]=func_Ab(z,I,rho,M,imask,S,II1,II2,JJ1,JJ2,Ivect)
    [nrows,ncols,nimgs] = size(I);
    % Depth gradient : 1st order forward diff
    zu = zeros(nrows,ncols);
    zv = zeros(nrows,ncols);
    zu(imask) = M(1:2:end-1,:)*z;
    zv(imask) = M(2:2:end,:)*z;

    % Norm of [-nabla z,1]
    dz2 = zu.^2+zv.^2+1;
    dz = sqrt(dz2);

    % A' field
    rho_over_d = transpose(rho(imask)./dz(imask));
    A1 = -bsxfun(@times,S(:,1),rho_over_d); % mxn
    A2 = -bsxfun(@times,S(:,2),rho_over_d); % mxn
    try
		A = sparse2([II1;II2],[JJ1(:);JJ2(:)],[A1(:);A2(:)]);
	catch
		A = sparse([II1;II2],[JJ1(:);JJ2(:)],[A1(:);A2(:)]);
	end
    % b' field
    b = Ivect(imask,:)' - bsxfun(@times,S(:,3),transpose(rho(imask)./dz(imask)));
    b = b(:);

func_f = 0.5*sum((A*M*z-b).^2)/(nimgs);
end



















