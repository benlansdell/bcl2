function [z_1, z_2, z_3] = solve_cubic(a_1, a_2, a_3, a_4)
	% solve_cubic		Find exact roots of cubic a_1*z^3 + a_2*z^2 + a_3*z + a_4 = 0
	%
	% Usage:
	%			[z_1, z_2, z_3] = cubic(a_1, a_2, a_3, a_4);
	%
	% Input:
	%			a_1 = z^3 coefficient
	%			a_2 = z^2 coefficient
	%			a_3 = z coefficient
	%			a_4 = constant
	%
	% Output:
	%			z_1 = First root
	%			z_2 = Second root (possibly repeated)
	%			z_3 = Third root (possibly repeated)
	% Still not working sadly...
	
	A = a_2/a_1;
	B = a_3/a_1;
	C = a_4/a_1;

	a = A/3;
	b = B/3;

	alpha = a^2 - b;
	beta = 2*a^3 - 3*a*b + C;

	if (alpha > 0)
		disp('alpha > 0')
		if (beta == 0)
			disp('beta == 0')
			z_1 = -a;
			z_2 = (3*alpha)^(0.5)-a;
			z_3 = -1*(3*alpha)^(0.5)-a;
		elseif (beta > 2*alpha^(3/2))
			disp('beta > 2*alpha^(3/2)');
			phi = acosh(abs(beta)/(2*(alpha^(3/2))))/3;
			z_1 = -2*alpha^(1/2)*cosh(phi)-a;
			z_2 = -sqrt(alpha)*cosh(phi)-a+i*sqrt(3*alpha)*sinh(phi);
			z_3 = -sqrt(alpha)*cosh(phi)-a+i*sqrt(3*alpha)*sinh(phi);
		elseif (beta < -2*alpha^(3/2))
			disp('beta < -2*alpha^(3/2)')
			phi = acosh(abs(beta)/(2*alpha^(3/2)))/3;			
			z_1 = -2*alpha^(0.5)*cosh(phi)-a;
			z_2 = -alpha^(0.5)*cosh(phi)-a-i*(3*alpha)^(0.5)*sinh(phi);
			z_3 = -alpha^(0.5)*cosh(phi)-a+i*(3*alpha)^(0.5)*sinh(phi);			
		else
			disp('|beta| < 2*alpha^(3/2)')
			phi = asin(beta/(2*alpha^(3/2)))/3;
			z_1 = 2*alpha^(0.5)*sin(phi)-a;
			z_2 = -2*alpha^(0.5)*sin(pi/3+phi)-a;
			z_3 = 2*alpha^(0.5)*sin(pi/3-phi)-a;
		end	
	elseif (alpha == 0)
		disp('alpha == 0')
		z_1 = -beta^(1/3)-a;
		z_2 = beta^(1/3)/2-a+3*i*beta^(2/3)/4;
		z_3 = beta^(1/3)/2-a-3*i*beta^(2/3)/4;
	elseif (alpha < 0)
		disp('alpha < 0')
		theta = asinh(beta/(2*(-alpha)^(3/2)))/3;
		z_1 = -2*(-alpha)^(0.5)*sinh(theta) - a;
		z_2 = (-alpha)^(0.5)*sinh(theta) - a + i*(-3*alpha)^(0.5)*cosh(theta);
		z_3 = (-alpha)^(0.5)*sinh(theta) - a - i*(-3*alpha)^(0.5)*cosh(theta);
	end
	
	%Check the answer is within absolute tolerance
	tol = 1e-6;
	f = @(z) a_1*z^3 + a_2*z^2 + a_3*z + a_4;
	if (max([abs(f(z_1)) abs(f(z_2)) abs(f(z_3))]) > tol)
		throw(MException('CubicError:OutsideTol', ['Computed roots outside forward error tolerence of ' num2str(tol)]));
	end
end