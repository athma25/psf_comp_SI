function output=nospace(id,run,rep,viz)

	tic
	close all


	% Checking output directory
	if ~exist(sprintf('output/%s',id));
		error('Batch folder does not exist');
	end

	%Reading parameters
	par=readmatrix(sprintf('output/%s/run%03d.txt',id,run));
	% Inputs: l - niche space [-l l]; m^2 - number of species; s0 - initial soil
	l=par(2,2);
	m=par(3,2);
	s0=par(4,2);

% Repetions with difference combinations
	rng(rep,'Threefry');

	% Parameters
	n=m^2;

	rp=NaN(1,n);		% Resource preference
	sp=NaN(n,1);		% Soil preference
	a=linspace(-l,l,m);
	dl=a(2)-a(1);
	for i=1:m
		rp(1,1+m*(i-1):m+m*(i-1))=a-0.5*dl+dl*rand(1,m);
		sp(1+m*(i-1):m+m*(i-1),1)=1.2*randn(m,1);		% Normally distributed
		%sp(1+m*(i-1):m+m*(i-1),1)=a(i)+0.5*rand(m,1);	% Uniformish random
	end

	%rp=linspace(-l,l,n);% Resource preference
	%sp=-l+2*l*rand(n,1);	% Soil preference
	%sp=-l/10+l*rand(n,1)/5;	% Soil preference
	%sp=zeros(n,1);	% Soil preference

%	rM=mvnrnd([0 0],[2 0; 0 1],n);
%	rp=rM(:,1)';
%	sp=rM(:,2);

	% Constant
	nu=par(5,2)*ones(n,1);				% Sensitivity to competition
	et=par(6,2)*ones(n,1);				% Conditioning strength
	sg_c=par(7,2)*ones(n,1);			% Range of competitive interactions
	sg_s=par(8,2)*ones(n,1);			% Range of selection

	% Bivariate log normal
%	m3=par(7,2);
%  m4=par(8,2);
%	v3=0.05;
%	v4=0.05;
%	c=0.5;
%	rM=exp(mvnrnd([log(m3^3/(m3^2+v3)) log(m4^3/(m4^2+v4))],[log(1+v3/m3^2) log(1+(c*sqrt(v3*v4))/(m3*m4)); log(1+(c*sqrt(v3*v4))/(m3*m4)) log(1+v4/m4^2)],n));
%	sg_c=rM(:,1);			% Range of competitive interactions
%	sg_s=rM(:,2);			% Range of selection

	et0=par(9,2);			% Strength of exogeneous driver

	imm=par(10,2);		% Immigration

	% Calculating competition coefficient alp
	rx=-20:0.01:20;
	
	% Linear infinite niche
 	ux=exp(-((rp-rx').^2)./(sg_c.^2)');		% Effective (resource)  utilization curve - dim(x) x n species
%	vx=exp(-((rp-rx').^4)./(81*sg_c.^4)');		% Total (resource)  utilization curve - dim(x) x n species

	% Circular niche with period 2l
%	ux=exp(-((rp-rx').^2)./(sg_c.^2)')+exp(-(((rp-2*l)-rx').^2)./(sg_c.^2)');		% Effective (resource)  utilization curve - dim(x) x n species
	vx=exp(-((rp-rx').^4)./(81*sg_c.^4)')+exp(-(((rp-2*l)-rx').^4)./(81*sg_c.^4)')+exp(-(((rp+2*l)-rx').^4)./(81*sg_c.^4)');		% Total (resource)  utilization curve - dim(x) x n species	

	alp=(ux'*vx)-ux(1,:)'*vx(1,:)/2-ux(end,:)'*vx(end,:)/2;
	alp_d=(ux'*ux)-ux(1,:)'*ux(1,:)/2-ux(end,:)'*ux(end,:)/2;
	alp=alp./diag(alp_d);

	% Analytic check when resource utilization ux is Gaussian
%	alp(end-9:end,end-9:end)
%	alp_est=exp(-(rp-rp').^2/(4*sg_c(1)^2));
%	alp_est(end-9:end,end-9:end)
%	find(alp(end-9:end,end-9:end)~=alp_est(end-9:end,end-9:end))

%	alp=exp(-(rp-rp').^4./(sg_c.^4));

% Numerics
	t0=par(12,2);
	tf=par(13,2);
	y0=[0.01*ones(n,1); s0];
	[t,y]=ode45(@(t,y) d_eq(n,y(1:n),y(n+1),sp,sg_s,alp,et,nu,et0,s0,imm),[t0 tf],y0);

	pop=y(end,1:n);
	s_end=y(end,end);

	% Figures
	fig1=figure('Visible',viz);
	scatter(rp,sp,100*ones(n,1),pop,'filled');
	xlim([-1.2*l 1.2*l]);
	ylim([-1.2*l 1.2*l]);
	cbar=colorbar;
	cbar.Label.String='Abundance';
	hold 'on';
	plot([-1.2*l 1.2*l],[s_end s_end],'-k');
	hold 'off';
	xlabel('Resource preference \zeta','FontSize',15);
	ylabel('Soil preference \epsilon','FontSize',15);
	title(sprintf('s(0)=%1.2g \t \\iota=%1.2g \t \\sigma_c=%1.2g \t \\sigma_s=%1.2g',s0,imm,par(7,2),par(8,2)),'FontWeight','normal');
	 
	fig2=figure('Visible',viz);
	scatter(rp,sp,100*ones(n,1),sg_c,'filled');
	xlim([-1.2*l 1.2*l]);
	ylim([-1.2*l 1.2*l]);
	cbar=colorbar;
	cbar.Label.String='Resource niche width';
	hold 'on';
	plot([-1.2*l 1.2*l],[s_end s_end],'-k');
	hold 'off';
	xlabel('Resource preference \zeta','FontSize',15);
	ylabel('Soil preference \epsilon','FontSize',15);
	title(sprintf('s(0)=%1.2g \t \\iota=%1.2g \t \\sigma_c=%1.2g \t \\sigma_s=%1.2g',s0,imm,par(7,2),par(8,2)),'FontWeight','normal');
	
	fig3=figure('Visible',viz);
	scatter(rp,sp,100*ones(n,1),sg_s,'filled');
	xlim([-1.2*l 1.2*l]);
	ylim([-1.2*l 1.2*l]);
	cbar=colorbar;
	cbar.Label.String='Soil niche width';
	hold 'on';
	plot([-1.2*l 1.2*l],[s_end s_end],'-k');
	hold 'off';
	xlabel('Resource preference \zeta','FontSize',15);
	ylabel('Soil preference \epsilon','FontSize',15);
	title(sprintf('s(0)=%1.2g \t \\iota=%1.2g \t \\sigma_c=%1.2g \t \\sigma_s=%1.2g',s0,imm,par(7,2),par(8,2)),'FontWeight','normal');

	fig4=figure('Visible',viz);
	stem(rp,pop);
	xlim([-1.2*l 1.2*l]);
	ylim([0 1]);
	xlabel('Resource preference \zeta');
	ylabel('Abundance');
	title(sprintf('s(0)=%1.2g \t \\iota=%1.2g \t \\sigma_c=%1.2g \t \\sigma_s=%1.2g',s0,imm,par(7,2),par(8,2)),'FontWeight','normal');

	fig5=figure('Visible',viz);
	plot(t,y(:,1:n),'-','LineWidth',1);
	ylim([0 1]);
	yyaxis right
	yaxr=gca;
	yaxr.YColor=[0 0 0];
	plot(t,y(:,end),'-k','LineWidth',2);
	ylim([-1.2*l 1.2*l]);
	xlabel('Time');
	ylabel('Density');
	title(sprintf('s(0)=%1.2g \t \\iota=%1.2g \t \\sigma_c=%1.2g \t \\sigma_s=%1.2g',s0,imm,par(7,2),par(8,2)),'FontWeight','normal');

	figure('Visible',viz)
	stem(sp,pop);
	xlim([-1.2*l 1.2*l]);
	xlabel('Soil preference \epsilon');
	ylabel('Abundance');

	writematrix([t,y],sprintf('output/%s/run%03d_%03d_out.txt',id,run,rep));
	output=[rp NaN NaN; sp' NaN NaN; t y];
	print(fig1,sprintf('output/%s/run%03d_%03d_2d.png',id,run,rep),'-dpng');
	print(fig2,sprintf('output/%s/run%03d_%03d_2d_c.png',id,run,rep),'-dpng');
	print(fig3,sprintf('output/%s/run%03d_%03d_2d_s.png',id,run,rep),'-dpng');
	print(fig4,sprintf('output/%s/run%03d_%03d_rp.png',id,run,rep),'-dpng');
	print(fig5,sprintf('output/%s/run%03d_%03d_dyn.png',id,run,rep),'-dpng');

	writematrix([s0 et0 imm l m tf s_end; rp' sp nu et sg_c sg_s pop'],sprintf('output/%s/run%03d_%03d_final.txt',id,run,rep));
	toc
	fprintf('Finished run%03d_%03d\n',run,rep);
end
