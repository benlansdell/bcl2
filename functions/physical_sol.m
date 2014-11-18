function sol = physical_sol(S, s)
	%For solution structure return only those which are 'physically' possible.
	%Those that are real and between 0 and 1...
	if (nargin < 2)
		s = 1
	end
	len = length(S.b);
	sol = [];
	for idx=1:len
		b = S.b(idx);
		t = S.t(idx);
		m = S.m(idx);
		d = S.d(idx);
		bh = S.bh(idx);
		th = S.th(idx);
		mh = S.mh(idx);
		b4 = S.b4(idx);
		if (isfield(S, 'x'))
			x = S.x(idx)
			xh = S.xh(idx);
			thx = S.thx(idx);
			if (physical(b, s) & physical(t, s) & physical(m, s) & physical(d, s) & physical(bh, s) & physical(th, s) & physical(mh, s) & physical(b4, s) & physical(x, s) & physical(xh, s) & physical(thx, s))
				sol = [sol; b t m d bh th mh b4 x thx xh];
			end	
		else
			if (physical(b, s) & physical(t, s) & physical(m, s) & physical(d, s) & physical(bh, s) & physical(th, s) & physical(mh, s) & physical(b4, s))
				sol = [sol; b t m d bh th mh b4];
			end
		end
	end

	size_sol = size(sol);
	if (size_sol(1) ~= 1)
		throw(MException('SolutionError:tooManySolutions', ['Too many/few physically possible solutions: ' num2str(size_sol(1)) ' instead of 1.']));
	end
	sol = double(sol);
end
