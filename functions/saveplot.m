function saveplot(h, filename)
%saveplot	Saves the specified plot
%
% Usage:
%			saveplot(h, filename)
%
% Input:
%			h = figure object. Can just use 'gcf'
%			filename = output filename. Write in encapsulated postscript (eps)
%
% Examples:
%			plot(0:10, (0:10).^2);
%			saveplot(gcf, './images/test.eps');
	set(h, 'paperunits', 'inches');
	set(h, 'papersize', [6 4]);
	set(h, 'paperposition', [0 0 6 4]);
	print(h,  filename, '-depsc');
end