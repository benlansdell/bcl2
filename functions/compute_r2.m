function [r2, bakmcl1r2, tbidmcl1r2] = compute_r2(observations, resnorm, f_tbid, time_points)
	%compute_r2		Function to compute coefficient of determination for
	%				fitted data.
	len = size(observations);
	if (len(1) == 13)

		s_tot = (observations - repmat(mean(observations,2),1,len(2))).^2;
		ss_tot = sum(sum(s_tot));
		r2 = 1-resnorm/ss_tot;

		%Compute bakmcl1 (tbid case) specific r2 value, row 8
		residuals = f_tbid(8,time_points) - observations(9,:);
		s_tot = (observations(9,:) - mean(observations(9,:))).^2;
		ss_tot = sum(s_tot);
		bakmcl1r2 = 1 - sum(residuals.^2)/ss_tot;

		%Compute tbidmcl1 (tbid case) specific r2 value, row 7
		residuals = f_tbid(7,time_points) - observations(13,:);
		s_tot = (observations(13,:) - mean(observations(13,:))).^2;
		ss_tot = sum(s_tot);
		tbidmcl1r2 = 1 - sum(residuals.^2)/ss_tot;

	else
		throw(MException('ValueError:invalidObservation', 'Invalid observation matrix: must be provided by one of observations_exp4_eqlm(), observations_exp4(), observations_exp4_multi()'));
	end
end