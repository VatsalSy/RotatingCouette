%% Getting the triple point contact line
% Author: Vatsal Sanjay
% vatsalsanjay@gmail.com
% Physics of Fluids
clc
clear
close
load('GerrisData.mat')
place = sprintf('lastNewt');
ll=evalc(sprintf('!./getUtheta %s', place));
lolo=textscan(ll,'%f %f\n');
rNewtBas=lolo{1}/0.25; uNewtBas=lolo{2}/0.25;
place = sprintf('lastShThn');
ll=evalc(sprintf('!./getUtheta %s', place));
lolo=textscan(ll,'%f %f\n');
rShThnBas=lolo{1}/0.25; uShThnBas=lolo{2}/0.25;
place = sprintf('lastShThn2');
ll=evalc(sprintf('!./getUtheta %s', place));
lolo=textscan(ll,'%f %f\n');
rShThnBas2=lolo{1}/0.25; uShThnBas2=lolo{2}/0.25;
% place = sprintf('lastHB');
% ll=evalc(sprintf('!./getUtheta %s', place));
% lolo=textscan(ll,'%f %f\n');
% rHBBas=lolo{1}/0.25; uHBBas=lolo{2}/0.25;
place = sprintf('lastBing');
ll=evalc(sprintf('!./getUtheta %s', place));
lolo=textscan(ll,'%f %f\n');
rBingBas=lolo{1}/0.25; uBingBas=lolo{2}/0.25;

eta = 2; % R2/R1 is eta
rNewtTH = linspace(1.0,eta,1001);
uNewtTh = (eta^2/(eta^2-1))*(1./rNewtTH - rNewtTH/eta^2);
rShThnTH = linspace(1.0,eta,1001);
uShThnTh = (eta^4/(eta^4-1))*(rShThnTH.^(-3) - rShThnTH/eta^4);
% rHBTH = linspace(0.25,0.5,1001);
% uHBTh = 0.25*(0.5*0.25/(0.5^2-0.25^2))*(0.5./rNewtTH - rNewtTH/0.5);
rBingTH = linspace(1.0,eta,1001);
tauy = 10/sqrt(2); V0 = 0.25; R1 = 0.25; mu0 = 1.0;
A = (tauy*R1)/(mu0*V0);
fun  = @(x) x^2 - log(x^2) - 1 - 2/A;
Ry = fzero(fun, 1.5);
sgn = -1;
uBingTh = rBingTH.*(1 + ((sgn*A*log(Ry) - 1)/(1 - Ry^(-2)))*(1 - rBingTH.^(-2)) - sgn*A*log(rBingTH));
for i = 1:1:length(rBingTH)
    if (rBingTH(i) > Ry)
        uBingTh(i) = 0.0;
    end
end
figure1 = figure('visible','on','WindowState','fullscreen','Color',[1 1 1]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');

plot(rNewtTH,uNewtTh,'r-','LineWidth',3,'DisplayName','Theory: Newtonian')
plot(rNewtGer/0.25,uNewtGer/0.25,'r.','MarkerSize',40,'DisplayName','Gerris: Newtonian')
plot(rNewtBas,uNewtBas,'r*','MarkerSize',15,'LineWidth',2,'DisplayName','Basilisk: Newtonian')

plot(rShThnTH,uShThnTh,'g-','LineWidth',3,'DisplayName','Theory: Power Law, n = 0.5')
plot(rShThnGer/0.25,uShThnGer/0.25,'g.','MarkerSize',40,'DisplayName','Gerris: Power Law, n = 0.5')
plot(rShThnBas,uShThnBas,'g*','MarkerSize',15,'LineWidth',2,'DisplayName','Basilisk: Power Law, n = 0.5')
plot(rShThnBas2,uShThnBas2,'k*','MarkerSize',15,'LineWidth',2,'DisplayName','Basilisk: Power Law, n = 0.38')
% plot(rHBGer,uHBGer,'c-','LineWidth',3,'DisplayName','Theory: Herschel-Bulkley')
% plot(rHBGer/0.25,uHBGer/0.25,'c.','MarkerSize',40,'DisplayName','Gerris: Herschel-Bulkley')
% plot(rHBBas,uHBBas,'c*','MarkerSize',15,'DisplayName','Basilisk: Herschel-Bulkley')

plot(rBingTH,uBingTh,'b-','LineWidth',3,'DisplayName','Theory: Bingham')
plot(rBingGer/0.25,uBingGer/0.25,'b.','MarkerSize',40,'DisplayName','Gerris: Bingham')
plot(rBingBas,uBingBas,'b*','MarkerSize',15,'LineWidth',2,'DisplayName','Basilisk: Bingham')

box(axes1,'on');
set(axes1,'DataAspectRatio',[1 1 1],'FontName','times new roman','FontSize',...
    30,'FontWeight','bold','LineWidth',3,'PlotBoxAspectRatio',[1 1 1]);
xlim([1.0 2.0])
ylim([0 1.0])
axis square
xlabel('\boldmath{$r/R_1$}','LineWidth',2,'FontWeight','bold','FontSize',50,...
            'FontName','times new roman',...
            'Interpreter','latex');
ylabel('\boldmath{$v_\theta$}','LineWidth',2,'FontWeight','bold','FontSize',75,...
    'FontName','times new roman',...
    'Interpreter','latex'); 

legend1 = legend(axes1,'show');
set(legend1,'LineWidth',3,'Interpreter','latex','FontSize',25,...
    'EdgeColor',[1 1 1]);
annotation(figure1,'textbox',...
    [0.495047619047627 0.458471760797336 0.271619047619044 0.0735807267180651],...
    'String','\boldmath{$\mu_{eq} = \frac{\tau_y}{2\|D_{ij}\|} + \mu_0\|D_{ij}\|^{n-1}$}',...
    'Interpreter','latex',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');
