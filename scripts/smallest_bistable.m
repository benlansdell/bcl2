%smallest_bistable.m
%Script to demonstrate the bistability phenomenona in 'smallest bistability' paper

%Set parameters
k1 = 8;
k2 = 1;
k3 = 1;
k4 = 1.5;
k5 = 0.6;
%k5 = 0;
S = 0.01:.01:2;

R = [];
for s = S
 R = [R, roots([-k3*k2/k1/s k2 -k4 k5])];
end
R

r1 = imag(R(1,:))==0;
r2 = imag(R(2,:))==0;
r3 = imag(R(3,:))==0;

%Plot exact and numerical solution to check they are the same...
plot(S(r1), R(1,r1), '-b', S(r2), R(2,r2), '--b', S(r3), R(3,r3), '-b');
xlim([0 2]);
ylim([0 15]);
xlabel('s');
ylabel('x*');

%Save image
saveplot(gcf, './images/demos/smallest_bistable.eps');
