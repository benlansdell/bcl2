function [b, t, m, d, c, bh, th, mh, ch, b4, x, thx, xh] = steadystate_conc(k, ic, p, u)
	%	
	%steadystate_conc	Function to compute the steady state concentrations for a
	%			set of kinetic parameters and initial conditions. Uses MATLAB's
	%			symbolic toolbox to compute analytically.
	%
	% Usage:
	%			[b, t, m, d, c, bh, th, mh, ch, b4] = steadystate_conc(k, ic)
	%			ss = steadystate_conc(k, ic);
	%
	% Input:
	%			k = 7x2 matrix of kinetic parameters. Can be taken and adjusted by
	%				load_model().
	%			ic = 10x1 vector of initial conditions
	%			p (optional) = 8x1 vector of protein production rates
	%			u (optional) = 8x1 vector of protein degradation rates
	%
	% Output:
	%			ss = 13x1 vector of the following species concentrations
	%			b = Bak steady state conc
	%			t = tBid/tBim steady state conc
	%			m = Mcl-1 steady state conc
	%			d = Bak* dimer steady state conc
	%			c = Cyt C mitos steady state conc
	%			bh = Bak* steady state conc
	%			th = tBid/tBim:Mcl-1 steady state conc
	%			mh = Bak*:Mcl-1 steady state conc
	%			ch = Cyt C S/N steady state conc
	%			b4 = Bak* multimer steady state conc
	%			x = Bcl-xl steady state conc
	%			thx = Bcl-xl:tBid/tBim steady state conc
	%			xh = Bcl-xl:Bak* steady state conc
	%
	% Examples:
	%			[k, ic] = load_model('tBim', 'reversible_exp4_BIAcore');
	%			[b, t] = steadystate_conc(k, ic);
	%			disp(['Steady-state tBim concentration is ' num2str(t) ' nM']);
	B0 = ic(1);
	T0 = ic(2);
	M0 = ic(3);
	C0 = ic(5);
	X0 = ic(13);

	%Provided that parameters are loaded from load_model() and adjusted by fit_exp4_*() then
	%we know the parameters are non-zero and that cytochrome c release is therefore
	%'inevitable'. Thus...
	c = 0;
	ch = C0;
	%Non-dimensionalise the system
	alpha = T0/M0;
	lambda = M0/B0;
	nu = X0/B0;
	kappa = T0/X0;
	
	K1D = k(1,2)/k(1,1)/T0;
	K2A = k(2,1)*B0/k(1,1)/T0;
	K2D = k(2,2)/k(1,1)/T0;
	K3A = k(3,1)*M0/k(1,1)/T0;
	K3D = k(3,2)/k(1,1)/T0;
	K4A = 2*k(4,1)*B0/k(1,1)/T0;
	K4D = k(4,2)/k(1,1)/T0;
	K5A = k(5,1)*B0/k(1,1)/T0;
	K6A = k(6,1)*B0/k(1,1)/T0;
	K6D = k(6,2)/k(1,1)/T0;
	K7A = k(7,1)*B0/k(1,1)/T0/4;
	K8A = k(8,1)*B0/k(1,1)/T0;
	K8D = k(8,2)/k(1,1)/T0;
	K9A = k(9,1)*X0/k(1,1)/T0;
	K9D = k(9,2)/k(1,1)/T0;

	%Solve dimensionless system first...very slow unfortunately
	syms b t m d bh th mh b4 x thx xh;
	if (X0 > 0)
		S = solve(t+th+thx-1, m+mh+alpha*th-1, b+bh+d+b4+lambda*mh+nu*xh-1, x+xh+kappa*thx-1, K2A*bh*m-K2D*mh, K3A*t*m-K3D*th, K6A*d*d-K6D*b4, K4A*bh*bh-K4D*d, b*t-K1D*bh+K5A*b*bh, K8A*x*bh-K8D*xh, K9A*t*x-K9D*thx);
	else
		S = solve(t+th-1, m+mh+alpha*th-1, b+bh+d+b4+lambda*mh-1, K2A*bh*m-K2D*mh, K3A*t*m-K3D*th, K6A*d*d-K6D*b4, K4A*bh*bh-K4D*d, b*t-K1D*bh+K5A*b*bh);	
	end
	sol = physical_sol(S);
	%Solve the dimensional system...
	%Solve it numerically?? I suppose...
	%if (X0 > 0)
	%	% 1  2  3  4  5  6   7   8   9   10  11 12   13
	%	%[b, t, m, d, c, bh, th, mh, ch, b4, x, thx, xh];
	%	options = optimset('TolX', 1e-12, 'TolFun', 1e-12, 'Display', 'iter');
	%	[sol, fval, exitflag] = fsolve(@(x)[x(2)+x(7)+x(12)-T0, x(3)+x(8)+x(7)-M0, x(1)+x(6)+2*x(4)+4*x(10)+x(8)+x(13)-B0, x(11)+x(13)+x(12)-X0, k(2,1)*x(6)*x(3)-k(2,2)*x(8), k(3,1)*x(2)*x(3)-k(3,2)*x(7), k(6,1)*x(4)*x(4)-k(6,2)*x(10), k(4,1)*x(6)*x(6)-k(4,2)*x(4), k(1,1)*x(1)*x(2)-k(1,2)*x(6)+k(5,1)*x(1)*x(6), k(8,1)*x(11)*x(6)-k(8,2)*x(13), k(9,1)*x(2)*x(11)-k(9,2)*x(12), x(5)+x(9)-C0], ic, options);
	%	%S = solve(t+th+thx-T0, m+mh+th-M0, b+bh+2*d+4*b4+mh+xh-B0, x+xh+thx-X0, k(2,1)*bh*m-k(2,2)*mh, k(3,1)*t*m-k(3,2)*th, k(6,1)*d*d-k(6,2)*b4, k(4,1)*bh*bh-k(4,2)*d, k(1,1)*b*t-k(1,2)*bh+k(5,1)*b*bh, k(8,1)*x*bh-k(8,2)*xh, k(9,1)*t*x-k(9,2)*thx);
	%else
	%	%S = solve(t+th-T0, m+mh+th-M0, b+bh+2*d+4*b4+mh-B0, k(2,1)*bh*m-k(2,2)*mh, k(3,1)*t*m-k(3,2)*th, k(6,1)*d*d-k(6,2)*b4, k(4,1)*bh*bh-k(4,2)*d, k(1,1)*b*t-k(1,2)*bh+k(5,1)*b*bh);
	%	options = optimset('TolX', 1e-12, 'TolFun', 1e-12, 'Display', 'iter');
	%	[sol, fval, exitflag] = fsolve(@(x)[x(2)+x(7)-T0, x(3)+x(8)+x(7)-M0, x(1)+x(6)+2*x(4)+4*x(10)+x(8)-B0, k(2,1)*x(6)*x(3)-k(2,2)*x(8), k(3,1)*x(2)*x(3)-k(3,2)*x(7), k(6,1)*x(4)*x(4)-k(6,2)*x(10), k(4,1)*x(6)*x(6)-k(4,2)*x(4), k(1,1)*x(1)*x(2)-k(1,2)*x(6)+k(5,1)*x(1)*x(6), x(5)+x(9)-C0], ic, options);
	%	%[sol, fval, exitflag] = solve(t+th-T0, m+mh+th-M0, b+bh+2*d+4*b4+mh-B0, k(2,1)*bh*m-k(2,2)*mh, k(3,1)*t*m-k(3,2)*th, k(6,1)*d*d-k(6,2)*b4, k(4,1)*bh*bh-k(4,2)*d, k(1,1)*b*t-k(1,2)*bh+k(5,1)*b*bh);
	%end
	%sol = physical_sol(S, 0);
	b = sol(1)*B0;
	t = sol(2)*T0;
	m = sol(3)*M0;
	d = sol(4)*B0/2;
	bh = sol(5)*B0;
	th = sol(6)*T0;
	mh = sol(7)*M0;
	b4 = sol(8)*B0/4;
	if (X0 > 0)
		x = sol(9);%*X0;
		thx = sol(10);%*T0;
		xh = sol(11);%*X0;
	else
		x = 0;
		thx = 0;
		xh = 0;
	end
	if (nargout == 1)
		%names = {'Bak'; 'tBim'; 'Mcl1'; 'DimBak'; 'CytCm'; 'aBak'; 'tBimMcl1'; 'aBakMcl1'; 'CytCc'; 'MulBak'; 'Bcl-xl'; 'tBimBclxl'; 'aBakBclxl'};
		b = [b, t, m, d, c, bh, th, mh, ch, b4, x, thx, xh];
	end
	%fval
	%exitflag
end