function plot_FLtau(tau, L, f, F)
	% view resolution
	stepL = 10;
	steptau = 10;

	figure();
	tiledlayout(1, 2);

	nexttile();
	waterfall(tau(1:steptau:end), L(1:stepL:end), f(1:stepL:end, 1:steptau:end));
	title('Andamento di f(\tau, L)');         
	xlabel('\tau [h]');                       
	ylabel('L [km]');                     
	zlabel('f(\tau, L)');  
	
	nexttile();
	waterfall(tau(1:steptau:end), L(1:stepL:end), F(1:stepL:end, 1:steptau:end));
	title('Andamento di f(\tau, L)');         
	xlabel('\tau [h]');                       
	ylabel('L [km]');                     
	zlabel('F(\tau, L)');
end