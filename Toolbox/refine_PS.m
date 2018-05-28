function [z_curr,rho_curr] = refine_PS(I,S,mask,lambda,z0,rho0,M, ...
    tol_alternate,tol_iPiano,maxit_alternate,maxit_iPiano, ...
    max_backtracking,L_init,eta,mu,c,delta_init,grad_approx)
	
    if nargin<18
        % Use approximated gradient
        grad_approx = true; 
        if nargin<17
            % delta on 0th iteration, must be >=c, should be >c
            delta_init = 1;
            if nargin<16
                % parameter to bound alpha away from the supremum
                c = 0.01;
                if nargin<15
                    % Denominator to decrease the initioal guess for the Lipschitz constant
                    mu = 1.05; 
                    if nargin<14
                        % Multiplier for the Lipschitz constant
                        eta = 1.2; 
                        if nargin<13
                            % Lipschit constant initialization
                            L_init = 1;
                            if nargin<12
                                % Bound on how often L can be increased in
                                % one iPiano iteration
                                max_backtracking = 100;
                                if nargin<11
                                    % Max. number of iterations for iPiano within each global iteration
                                    maxit_iPiano = 100;
                                    if nargin<10
                                        % Max. number of global iterations
                                        maxit_alternate = 500;
                                        if nargin<9
                                            % Tolerance for the outer, alternating optimisation
                                            tol_iPiano = 1e-8;
                                            
                                            if nargin<8
                                                % Tolerance for the outer, alternating optimisation
                                                tol_alternate = 1e-8;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
	
	% Aux. variables
	imask = find(mask>0);
	npix = length(imask);
	rho_curr = rho0; % Current rho-estimate
	z_curr = z0; % Current z-estimate
	z_next = z_curr; % Next z-estimate
	z_prev = z_curr; % Previous z-estimate
	[nrows,ncols,nimgs] = size(I); % Size of I
	Ivect = reshape(I,nrows*ncols,nimgs);
	II1 = transpose(1:nimgs*npix);
	JJ1 = 1:2:2*npix-1;
	JJ1 = repmat(JJ1,[nimgs 1]);
	II2 = transpose(1:nimgs*npix);
	JJ2 = 2:2:2*npix;
	JJ2 = repmat(JJ2,[nimgs 1]);

	
	% Current gradient
	[nrows,ncols,~] = size(I);
	zu_curr =  zeros(nrows,ncols);
	zv_curr =  zeros(nrows,ncols);
	zu_curr(imask) = M(1:2:end-1,:)*z_curr(imask);
	zv_curr(imask) = M(2:2:end,:)*z_curr(imask);

	func_f_curr = func_f(z_curr,I,rho_curr,M,imask,S,II1,II2,JJ1,JJ2,Ivect,grad_approx);	
	energy_curr_alternate = func_f_curr+0.5*lambda*sum((z_curr(imask)-z0(imask)).^2);


	% Albedo update (linear least squares)
	rho_curr = update_rho(I,zu_curr,zv_curr,S,mask);
		
	% Alternate optimization steps
	for it_alternate = 1:maxit_alternate

		% Reset Lipschitz constant estimate
		L = L_init;
        
        % Reset initial delta for the computation of beta
        delta = delta_init;

		energy_curr = energy_curr_alternate;

		
		% Depth update (iPiano)
		for it_iPiano = 1:maxit_iPiano

			% Compute f and its gradient at current estimate
			[func_f_curr,~,~,grad_f_curr] = func_f(z_curr,I,rho_curr,M,imask,S,II1,II2,JJ1,JJ2,Ivect,grad_approx);
            % store the quality of the approximated gradient

			% Lazy backtracking to set stepsize
			lc = 0; % lc = 1 if Lipschitz constant L is big enough
            while(lc < max_backtracking)
                % save the actually used Lipschitz constant for later
                L_curr = L;
                % supplemental variable b
                b = (delta+0.5*L)/(c+0.5*L);
                % see proof of Lemma 4.6 in Ochs et al 2014, +c is added to
                % enforce (b-1)/(b-0.5)>beta and thereby a descent of delta
                beta = (b-1)/(b+c-0.5);
				alpha = (1-beta)/(c+0.5*L); % Descent stepsize
						
				z_next(imask) = prox_g(z_curr(imask)-alpha*grad_f_curr+beta*(z_curr(imask)-z_prev(imask)),alpha,lambda,z0(imask)); % Compute next extimate

				z_dist = z_next - z_curr; % Evaluate the difference between current and next estimate	
				[func_f_next,zu_next,zv_next] = func_f(z_next,I,rho_curr,M,imask,S,II1,II2,JJ1,JJ2,Ivect,grad_approx); % Evaluate at the next estimate
				
				% Lipschitz test
                if(func_f_next <= func_f_curr+grad_f_curr'*z_dist(imask)+0.5*L*z_dist(imask)'*z_dist(imask))
					L = L/mu; % Advised by Ochs et al for speedup
					lc = max_backtracking; % if Lipschitz => stepsize is small enough
				else
					lc = lc+1; % if not Lipschitz => try smaller stepsize
					L = eta*L;
                end
            end
            delta = 1/alpha-L_curr/2-beta/(2*alpha);
		
			% Update auxiliary variables
			z_prev = z_curr;
			z_curr = z_next;
			zu_curr = zu_next;
			zv_curr = zv_next;

			% Compute full energy
			energy_next = func_f_next+0.5*lambda*sum((z_curr(imask)-z0(imask)).^2);
            
            
            % CV test
            %%% relative_residual = norm(z_curr(imask)-z_prev(imask))./norm(z_curr(imask));
            relative_residual = abs(energy_next-energy_curr)/abs(energy_curr);
            %~ disp(sprintf('It. %d.%d - res : %.7f - E(z) : %.7f - L : %.2f  ',it_alternate,it_iPiano,relative_residual,energy_next,L));
            if(relative_residual < tol_iPiano)
                break;
            end
			
			% Store current estimate for CV evaluation
			energy_curr = energy_next;			
		end

		% Albedo update (linear least squares)
		rho_curr = update_rho(I,zu_curr,zv_curr,S,mask); 



		% Global CV test
		relative_residual_alternate = abs(energy_next-energy_curr_alternate)/abs(energy_curr_alternate);

		fprintf('iPiano loop  %4d/%d ended in iteration %3d/%d - residual : %.9f - E(z) : %10.7f - L : %.2e\n', ...
            it_alternate,maxit_alternate,it_iPiano,maxit_iPiano, ...
            relative_residual_alternate,energy_next,L);


		if(relative_residual_alternate < tol_alternate)
			break;
		end
		energy_curr_alternate = energy_next;		
    end

end
























