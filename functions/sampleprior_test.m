function rd = sampleprior_test(mu, s1)
	stddev = [s1, s1];
	nz = size(mu, 2);
	zs = randn(1, nz);
	rd = mu + zs.*stddev;
end