function lp = prior(k, mu, s1)
	z = [k(1), k(2)];
	stddev = [s1, s1];
	p = normpdf(z, mu, stddev);
	lp = sum(log(p));
end