function obs = observations_exp6_multi()
	% observations_exp6_multi 		Return matrix of densitometry data to fit experiment 6 model to. Each 
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
	%		obs = observations_exp6_multi()
	%
	% Output:
	%		obs = 14x15 matrix of concentration values for 14 sets of observations over 15 time points.
	%
	% Examples:
	% 		obs = observations_exp6_multi()
	
	obs = [ 0.043105305	0.000693856	0.09873329	0.01417368	0.282553204	0.999623494	1.622121962	2.516921961	2.714091981	3.270253467	3.414746138	3.554614254	3.586706398	4.143876929	4.508449773;	%tBIM: Cyt C S/N release
			4.956894695	4.999306144	4.90126671	4.98582632	4.717446796	4.000376506	3.377878038	2.483078039	2.285908019	1.729746533	1.585253862	1.445385746	1.413293602	0.856123071	0.491550227;	%tBIM: Cyt C mitos
			0.781	0.846	2.577	10.143	4.027	2.564	1.600	2.899	0.073	0.775	0.748	0.642	0.991	0.293	0.358;				%tBIM: Bak*:Mcl-1
			2			2			2			2			2			2			2			2			2	2	2	2	2	2	2;				%tBIM: Bak* free
			2			2			2			2			2			2			2			2			2	2	2	2	2	2	2;				%tBIM: Bak* dimer
			2			2			2			2			2			2			2			2			2	2	2	2	2	2	2;				%tBIM: Bak* multimer
            6.859273759	11.04637059	13.19871582	20			20			20			10.71630187	19.55938416	20	20	20	20	20	20	19.64216667;				%tBIM: tBIM:Mcl-1
            0.003840869	0.010521271	0.006294496	0.018543081	0.028324882	0.071314757	0.097875765	0.233443547	0.314826241	0.376005415	0.450318083	0.455888639	0.459954039	0.899681844	1.057177511;	%tBid: Cyt C S/N release
			4.996159131	4.989478729	4.993705504	4.981456919	4.971675118	4.928685243	4.902124235	4.766556453	4.685173759	4.623994585	4.549681917	4.544111361	4.540045961	4.100318156	3.942822489;	%tBid: Cyt C mitos
			0.649	0.764	0.823	3.424	5.260	4.368	4.191	6.917	6.452	6.439	6.871	5.558	2.588	5.847	6.769;	%tBid: Bak*:Mcl-1
			1			1			1			1			1			1			1			1			1	1	1	1	1	1	1;				%tBid: Bak* free
			1			1			1			1			1			1			1			1			1	1	1	1	1	1	1;				%tBid: Bak* dimer
			0.5			0.5			0.5			0.5			0.5			0.5			0.5			0.5			0.5 0.5	0.5	0.5	0.5	0.5	0.5;			%tBid: Bak* multimer -- same form as Cyt C S/N release but rescaled
			0.99561264	1.076035389	2.042208005	2.70293412	3.890831606	3.628270574	5.628992636	11.52109612	9.902594541	8.56284851	12.55884448	13.87699783	9.534669212	15.40522999	13.23116667];			%tBid: tBid:Mcl-1
end