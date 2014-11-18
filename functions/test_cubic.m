%Script to test each case of cubic.m solver.
%Function performs a self check each time so we only need to run the function under each case
%and ensure no exceptions are raised.

%Set a, b, C in each case

disp('alpha > 0; beta == 0 case...');
a = 4;
b = 7;
C = 3*a*b - 2*a^3;
solve_cubic(1, 3*a, 3*b, C);
disp('alpha > 0; beta == 0 case success');

disp('alpha > 0; |beta| < 2*alpha^(3/2) case 1...');
a = 4;
b = 7;
C = 5;
solve_cubic(1, 3*a, 3*b, C);
disp('alpha > 0; |beta| < 2*alpha^(3/2) case 1 success');

disp('alpha > 0; |beta| < 2*alpha^(3/2) case 2...');
a = 4;
b = 7;
C = -88;
solve_cubic(1, 3*a, 3*b, C);
disp('alpha > 0; |beta| < 2*alpha^(3/2) case 2 success');

disp('alpha > 0; beta > 2*alpha^(3/2) case...');
a = 4;
b = 7 ;
C = 12;
solve_cubic(1, 3*a, 3*b, C);
disp('alpha > 0; beta > 2*alpha^(3/2) case success');

%disp('alpha > 0; beta < -2*alpha(3/2) case...');
%a = 4;
%b = 7;
%C = -100;
%solve_cubic(1, 3*a, 3*b, C);
%disp('alpha > 0; beta < -2*alpha(3/2) case success');

%disp('alpha == 0 case...');
%a = 4;
%b = 16;
%C = 20;
%solve_cubic(1, 3*a, 3*b, C);
%disp('alpha == 0 case success');

disp('alpha < 0 case...');
a = 4;
b = 32;
C = 20;
solve_cubic(1, 3*a, 3*b, C);
disp('alpha < 0 case success');

disp('Test successful.');