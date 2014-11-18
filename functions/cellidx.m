function idx = cellidx(cellArray, string)
	% cellidx		Search given cell array for string. Return linear index in array if found
	%
	% Usage:
	%			idx = cellidx(cellArray, string)
	%
	% Input:
	%			cellArray = cell array in which to search
	%			string = string to search for
	%
	% Output:
	%			idx = index within array. If string not found, return 0
	%
	% Example(s):
	%			j = cellidx({'1', '3', '2'}, '3');
	
	idx = find(ismember(cellArray, string)==1);
	if (length(idx) == 0)
		idx = 0;
	end
end