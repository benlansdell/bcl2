function [obs, obs_std] = observations_exp46()
	% observations_exp46 		Return matrix of densitometry data to fit experiment 14 model to. Each 
	%							row is a species, each column is a time point. This data set is obtained
	%							through a mix of rough estimates based on densitometry data and some 
	%							qualitative estimates of what should be happening in the MLM system -- 
	%							i.e. later Bak dimerisation and hence cytochrome c release. The tBid
	%							densitometry data is scaled assuming that at its peak almost all Bak has
	%							has been activated and is bound to Mcl-1 (assumed 9nM). By the end (180min)
	%							it is assumed almost all the Mcl-1 is bound to something and therefore the 
	%							tBid:Mcl-1 conc can be inferred from the Bak:Mcl-1 conc at that time (13nM).
	%
	% Usage:
	%		obs = observations_exp46()
	%
	% Output:
	%		obs = 14x12 matrix of concentration values for 14 sets of observations over 12 time points.
	%
	% Examples:
	% 		obs = observations_exp46()

	%Without pseudoobservations (as percentages, needs to be scaled)
	obs = [ ...
		6.3, 2.6, 4.5, 6.5, 15.3, 20.5, 26.9, 37.4, 64.2, 96.0; 							%BidBim: CytC 		(CytC limited. x 5/100)
		9.9, 8.4, 11.8, 36.9, 95.6, 97.1, 95.2, 72.6, 65.9, 26.1; 							%Bak:Mcl1 			(Bak limited. x B_0/100)
		100.0, 92.0, 83.3, 48.2, 5.0, 2.3, 2.4, 2.1, 2.7, 2.7; 								%Inactive Bak 		(Bak limited. x B_0/100)
		0.03*100/35 0.27*100/35 27*100/35	100	  100	100    100 	100	  100	100; 		%tBim:Mcl-1 		(Mcl-1 limited. x 35/100)
		1.0, 4.4, 5.6, 11.7, 14.1, 22.9, 24.0, 23.3, 21.9, 31.6; 							%Bid: CytC release	(CytC limited. x 5/100)
		11.8, 9.9, 19.5, 78.2, 91.8, 78.6, 86.1, 92.1, 67.4, 63.9; 							%Bak:Mcl1			(Bak limited. x B_0/100)
		100.0, 94.9, 64.1, 11.1, 2.9, 1.5, 1.2, 1.1, 1.2, 1.0;								%Inactive Bak		(Bak limited. x B_0/100)
		0.03*100/35 0.27*100/35 27*100/35	100	  100	100    100 	100	  100	100]; 		%tBid:Mcl-1			(Mcl-1 limited. x 35/100)

	%With pseudo-observations designed to keep activated bak complexes from forming
	%in Bid case and to make sure they form in BidBim case... try the model without this...
	%	obs = [ ...
	%6.3, 2.6, 4.5, 6.5, 15.3, 20.5, 26.9, 37.4, 64.2, 96.0; 							%BidBim: CytC 		(CytC limited. x 5/100)
	%9.9, 8.4, 11.8, 36.9, 95.6, 97.1, 95.2, 72.6, 65.9, 26.1; 							%Bak:Mcl1 			(Bak limited. x B_0/100)
	%2	 2	  2		2	  2	    2	  2	    2	  2	    2; 								%Bak* 				(Bak limited. x B_0/100)
	%2	 2	  2		2	  2	    2	  2	    2	  2	    2; 								%Bak*^2 			(Bak limited. x B_0/100)
	%2	 2	  2		2	  2	    2	  2	    2	  2	    2; 								%Multi Bak* 		(Bak limited. x B_0/100)
	%0.03*100/35 0.27*100/35 27*100/35	100	  100	100    100 	100	  100	100; 		%tBim:Mcl-1 		(Mcl-1 limited. x 35/100)
	%1.0, 4.4, 5.6, 11.7, 14.1, 22.9, 24.0, 23.3, 21.9, 31.6; 							%Bid: CytC release	(CytC limited. x 5/100)
	%11.8, 9.9, 19.5, 78.2, 91.8, 78.6, 86.1, 92.1, 67.4, 63.9; 							%Bak:Mcl1			(Bak limited. x B_0/100)
	%0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5;											%Bak*				(Bak limited. x B_0/100)
	%0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5;											%Bak^*				(Bak limited. x B_0/100)
	%0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5; 											%Multi Bak*			(Bak limited. x B_0/100)
	%0.03*100/35 0.27*100/35 27*100/35	100	  100	100    100 	100	  100	100]; 		%tBid:Mcl-1			(Mcl-1 limited. x 35/100)

	%Stdevs (as percentages, needs to be scaled)
	obs_std = [ ...
		4.7, 3.3, 3.1, 8.9, 3.3, 2.6, 5.9, 8.5, 11.6, 3.1; 				%%BidBim:CytC
		4.4, 4.1, 6.7, 3.0, 4.1, 2.1, 5.9, 4.6, 2.4, 5.7; 				%Bak:Mcl-1
		0.1, 6.4, 19.2, 21.7, 0.4, 0.7, 1.0, 0.6, 0.9, 0.5;				%Bak
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5; 									%BidBim:Mcl-1
		1.3, 5.0, 2.6, 5.9, 13.0, 9.1, 4.5, 4.4, 0.9, 0.1; 				%Bid: CytC
		3.4, 8.4, 11.4, 6.1, 11.6, 18.0, 9.7, 11.2, 7.8, 2.7; 			%Bak:Mcl1
		0.1, 0.2, 4.5, 0.2, 0.5, 0.7, 0.7, 0.5, 0.6, 0.8; 				%Bak
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5];									%Bid:Mcl-1

end

%Bid					
%cytochrome c release		Bak activation (PK)		Bak:Mcl-1 (IP)	
%Mean						Mean					Mean	
%1.0						0.0						11.8	
%4.4						5.1						9.9		
%5.6						35.9					19.5	
%11.7						88.9					78.2	
%14.1						97.1					91.8	
%22.9						98.5					78.6	
%24.0						98.8					86.1	
%23.3						98.9					92.1	
%21.9						98.8					67.4	
%31.6						99.0					63.9	
%
%BidBim					
%cytochrome c release		Bak activation (PK)		Bak:Mcl-1 (IP)	
%Mean						Mean					Mean
%6.3						0.0						9.9	
%2.6						8.0						8.4	
%4.5						16.7					11.8
%6.5						51.8					36.9
%15.3						95.0					95.6
%20.5						97.7					97.1
%26.9						97.6					95.2
%37.4						97.9					72.6
%64.2						97.3					65.9
%96.0						97.3					26.1
%
%Bid					
%cytochrome c release		Bak activation (PK)		Bak:Mcl-1 (IP)	
%Std						Std						Std
%1.3						0.0						3.4
%5.0						0.2						8.4
%2.6						4.5						11.4
%5.9						0.2						6.1
%13.0						0.5						11.6
%9.1						0.7						18.0
%4.5						0.7						9.7
%4.4						0.5						11.2
%0.9						0.6						7.8
%0.0						0.8						2.7
%
%BidBim					
%cytochrome c release		Bak activation (PK)		Bak:Mcl-1 (IP)	
%Std						Std						Std
%4.7						0.0						4.4
%3.3						6.4						4.1
%3.1						19.2					6.7
%8.9						21.7					3.0
%3.3						0.4						4.1
%2.6						0.7						2.1
%5.9						1.0						5.9
%8.5						0.6						4.6
%11.6						0.9						2.4
%3.1						0.5						5.7
