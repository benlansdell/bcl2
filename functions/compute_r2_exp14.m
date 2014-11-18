function [r2, bakmcl1r2bid, tbidmcl1r2, bakmcl1r2bim, tbimmcl1r2] = compute_r2_exp14(observations, resnorm, f_tbid, f_tbim, time_points)
	%compute_r2		Function to compute coefficient of determination for
	%				fitted data.
	len = size(observations);
	if (len(1) == 14)

		s_tot = (observations - repmat(mean(observations,2),1,len(2))).^2;
		ss_tot = sum(sum(s_tot));
		r2 = 1-resnorm/ss_tot;

		%Compute bakmcl1 (tbid case) specific r2 value, row 8
		residuals = f_tbid(8,time_points) - observations(10,:);
		s_tot = (observations(10,:) - mean(observations(10,:))).^2;
		ss_tot = sum(s_tot);
		bakmcl1r2bid = 1 - sum(residuals.^2)/ss_tot;

		%Compute tbidmcl1 (tbid case) specific r2 value, row 7
		residuals = f_tbid(7,time_points) - observations(14,:);
		s_tot = (observations(14,:) - mean(observations(14,:))).^2;
		ss_tot = sum(s_tot);
		tbidmcl1r2 = 1 - sum(residuals.^2)/ss_tot;

		%Compute bakmcl1 (tbim case) specific r2 value, row 8
		residuals = f_tbim(8,time_points) - observations(3,:);
		s_tot = (observations(3,:) - mean(observations(3,:))).^2;
		ss_tot = sum(s_tot);
		bakmcl1r2bim = 1 - sum(residuals.^2)/ss_tot;

		%Compute tbimmcl1 (tbim case) specific r2 value, row 7
		residuals = f_tbim(7,time_points) - observations(7,:);
		s_tot = (observations(7,:) - mean(observations(7,:))).^2;
		ss_tot = sum(s_tot);
		tbimmcl1r2 = 1 - sum(residuals.^2)/ss_tot;

	else
		throw(MException('ValueError:invalidObservation', 'Invalid observation matrix: must be provided by one of observations_exp14_multi()'));
	end
end
