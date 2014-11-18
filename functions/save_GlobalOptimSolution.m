function save_GlobalOptimSolution(globalsols, fn_out)
	%compute_r2		Function to compute coefficient of determination for
	%				fitted data.

	nrows = length(globalsols);
	ncols = 1+length(globalsols(1).X);
	csv = zeros(nrows, ncols);
	for idx = 1:nrows
		csv(idx,1) = globalsols(idx).Fval;
		csv(idx,2:end) = globalsols(idx).X;
	end
	%Write to file
	xlswrite(fn_out, csv);
end