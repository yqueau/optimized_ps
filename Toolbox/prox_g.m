function p = prox_g(z,alpha,lambda,z0)

	p = (z+alpha*lambda*z0)./(1+alpha*lambda);
end
