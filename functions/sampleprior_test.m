function rd = sampleprior_test(mu, s1)
	rd = zeros(size(mu));
	neg = 1:length(mu);
	while length(neg) > 0
		stddev = s1*ones(size(neg));
		nz = length(neg);
		zs = randn(1, nz);
		rd(neg) = mu(neg) + zs.*stddev;
		neg = find(rd<0);
	end
end