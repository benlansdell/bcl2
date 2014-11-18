function obs = observations_exp4_eqlmest()
	% observations_exp4_equlmest 		Return matrix of densitometry data to fit experiment 4 model to. Each 
	%							row is a species, each column is a time point. This data set is obtained
	%							through a mix of rough estimates based on densitometry data and some 
	%							qualitative estimates of what should be happeningin the MLM system -- 
	%							i.e. later Bak dimerisation and hence cytochrome c release. The
	%							densitometry data is scaled according to estimates of 'equilibruim'
	%							concentrations being reached after 90 minutes, based on initial
	%							concentrations and Biacore data. Since Biacore was ultimately decided
	%							to be unreliable this is a less than ideal solution.
	%
	% Usage:
	%		obs = observations_exp4_equlmest()
	%
	% Output:
	%		obs = 13x9 matrix of concentration values for 13 sets of observations over 9 time points.
	%
	% Examples:
	% 		obs = observations_exp4_eqlmest()
	
	obs = [0.152551289	0.131509732	0.135642895	0.164950778	0.166829488	0.19050124	0.192004208	0.225069512	0.680093184;	%tBIM: Cyt C S/N release
			3.794995115	4.321034042	4.321034042	4.321034042	3.832569324	4.80949876	4.734350342	4.696776133	2.780491471;	%tBIM: Cyt C mitos
			0		0		0		0		0		0		0		0		0; 												%tBIM: Bak*:Mcl-1
			0		0		0		0		0		0		0		0		0;												%tBIM: Bak* free
			0.224309392	0.193370166	0.199447514	0.242541436	0.245303867	0.280110497	0.282320442	0.330939227	1;				%tBIM: Bak* dimer -- same form as Cyt C S/N release but rescaled
			5		5		5		5		5		5		5		5		5;												%tBIM: later Bak* dimer -- 120 minutes as after -- force cytochrome c release
			0.132578077	0.152173913	0.16442131	0.184935701	0.253521127	0.330679731	0.282914881	0.345988977	0.407225964;	%tBid: Cyt C S/N
			3.123086344	3.919167177	4.378444581	3.85793019	3.919167177	4.500918555	4.654011023	4.654011023	3.704837722;	%tBid: Cyt C mitos
			0.708039462	0.735271749	1.2050287	4.186964126	6.944233184	6.358739013	6.944233184	5.534962332	4.5546; 		%tBid: Bak*:Mcl-1
			0		0		0		0		0		0		0		0		0;												%tBid: Bak* free
			0.32556391	0.373684211	0.403759398	0.454135338	0.622556391	0.812030075	0.694736842	0.84962406	1;				%tBid: Bak* dimer -- same form as Cyt C S/N release but rescaled
			1		1		1		1		1		1		1		1		1;												%tBid: later Bak* dimer -- 120 minutes as after -- force no cytochrome c release
			0.389981643	0.361489833	0.688255274	0.988309643	1.745123333	2.270441071	7.318833571	18.07449167	14.9582]; 		%tBid: tBid:Mcl-1
end