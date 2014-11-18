function obs = observations_exp4()
	% observations_exp4 		Return matrix of densitometry data to fit experiment 4 model to. Each 
	%							row is a species, each column is a time point. This data set is obtained
	%							through a mix of rough estimates based on densitometry data and some 
	%							qualitative estimates of what should be happeningin the MLM system -- 
	%							i.e. later Bak dimerisation and hence cytochrome c release. The tBid
	%							densitometry data is scaled assuming that at its peak almost all Bak has
	%							has been activated and is bound to Mcl-1 (assumed 9nM). By the end (80min)
	%							it is assumed almost all the Mcl-1 is bound to something and therefore the 
	%							tBid:Mcl-1 conc can be inferred from the Bak:Mcl-1 conc at that time (13nM).
	%
	% Usage:
	%		obs = observations_exp4()
	%
	% Output:
	%		obs = 13x9 matrix of concentration values for 13 sets of observations over 9 time points.
	%
	% Examples:
	% 		obs = observations_exp4()
	
	obs = [ 0.152551289	0.131509732	0.135642895	0.164950778	0.166829488	0.19050124	0.192004208	0.225069512	0.680093184;	%tBIM: Cyt C S/N release
			3.794995115	4.321034042	4.321034042	4.321034042	3.832569324	4.80949876	4.734350342	4.696776133	2.780491471;	%tBIM: Cyt C mitos
			0			0			0			0			0			0			0			0			0;				%tBIM: Bak*:Mcl-1
			0			0			0			0			0			0			0			0			0;				%tBIM: Bak* free
			0.224309392	0.193370166	0.199447514	0.242541436	0.245303867	0.280110497	0.282320442	0.330939227	1; 				%tBIM: Bak* dimer
			5			5			5			5			5			5			5			5			5;				%tBIM: later Bak* dimer -- 120 minutes as after -- force cytochrome c release
			0.132578077	0.152173913	0.16442131	0.184935701	0.253521127	0.330679731	0.282914881	0.345988977	0.407225964;	%tBid: Cyt C S/N release
			3.123086344	3.919167177	4.378444581	3.85793019	3.919167177	4.500918555	4.654011023	4.654011023	3.704837722;	%tBid: Cyt C mitos
			0.917647059	0.952941176	1.561764706	5.426470588	9			8.241176471	9			7.173529412	5.902941176;	%tBid: Bak*:Mcl-1
			0			0			0			0			0			0			0			0			0;				%tBid: Bak* free
			0.32556391	0.373684211	0.403759398	0.454135338	0.622556391	0.812030075	0.694736842	0.84962406	1;				%tBid: Bak* dimer -- same form as Cyt C S/N release but rescaled
			1			1			1			1			1			1			1			1			1;				%tBid: later Bak* dimer -- 120 minutes as after -- force no cytochrome c release
			0.338928571	0.314166667	0.598154762	0.858928571	1.516666667	1.973214286	6.360714286	15.70833333	13];			%tBid: tBid:Mcl-1
end