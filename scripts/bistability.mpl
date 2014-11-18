#Script to find steady-states of simplified ODE system only containing tBid, Bak, Bak* and Mcl-1. 
#No dimers or multimers...
with(plots);
T0 := 30;
M0 := 20;
B0 := 10;
Th := k5*T*M;
Mh := k4*Bh*M;
T := T0/(1+k5*M);
B := B0 - Bh - Mh;
Bh_exp := Bh = k1*B*Bh + k2*T*B;
B := B0 - Bh - Mh;
Bh_sols := {solve(Bh_exp, Bh)};
Bh := Bh_sols[2];
M_exp := M0 = M + Mh + Th;
M2 := M_exp - M - T0*k5*M/(1+k5*M);
M3 := M2*2*k1;
M4 := M3*(k4*M^2*k5+1+k5*M + k4*M);
M5 := M4/M/k4;
M6 := -B0*k1+M5+k5*M-10*k1*k5*M+1+T0*k2+T0*k2*k4*M;
M7 := lhs(M6)^2 = rhs(M6)^2;
M8 := simplify(M7);
M9 := M8*k4*k4*M*M;
M10 := collect(M9,M);
M11 := isolate(M10, M^6); 
M12 := collect(M11, M);
k1:=1;
k2:=2;
k4:=3;
k5:=4;
plotsetup(cps, plotoutput=cat(`./images/test.ps`), plotoptions=`noborder,width=8in,height=6in`);
print(display(plot(lhs(M12)-rhs(M12), M=-1..21, gridlines=true, labels=['T0', 'd']),symbolsize=1,axes=boxed));
#I don't think it's worth continuing...