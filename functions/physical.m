function is_phys = physical(sym_a, scaled)
	if (nargin < 2)
		scaled = 1
	end
	a = double(sym_a);
	if scaled
		if (imag(a) == 0 & real(a) >= 0 & real(a) <= 1)
			is_phys = 1;
		else
			is_phys = 0;
		end
	else
		if (imag(a) == 0 & real(a) >= 0)
			is_phys = 1;
		else
			is_phys = 0;
		end	
	end
end
