function ssr = ll_test(k, t, x0, options, data, sig, solver, pts)
	dfk = @(t,x) [k(2)*x(3)-k(1)*x(1)*x(2)*x(2); 2*k(2)*x(3)-2*k(1)*x(1)*x(2)*x(2); k(1)*x(1)*x(2)*x(2) - k(2)*x(3)];
	[T, f] = solver(dfk, t, x0, options);
	res = f(pts,:) - data;
	ssr = -sum(sum(res.^2./sig.^2))/2;
end