function lp = prior(k, mu, s1)
	if any(k< 0)
		lp = -inf;
	else
		z = [k(1), k(2)];
		stddev = [s1, s1];
		p = normpdf(z, mu, stddev);
		lp = sum(log(p));
	end
end