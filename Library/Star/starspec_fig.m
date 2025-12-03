% make figure for HZ

% sun
x=importdata('solar.txt');

% Create figure
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');

semilogy(x(:,1),x(:,2),'LineWidth',2,'DisplayName','Sun');

% Create ylabel
ylabel('Flux [W/m^2/nm]');
xlabel('Wavelength [nm]');
xlim(axes1,[100 500]);
ylim(axes1,[1e-5,1e+1]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',20,'XMinorTick','on','YMinorTick','on','YScale','log');

% bolometric flux
ls = trapz(x(:,1),x(:,2)); % w/m2

% other stars

x=importdata('hd97658.txt');
l=trapz(x(:,1),x(:,2));
x(:,2)=x(:,2)*ls/l;
semilogy(x(:,1),x(:,2),'LineWidth',2,'DisplayName','HD 97658 (K1V)');

x=importdata('eps_eri.txt');
l=trapz(x(:,1),x(:,2));
x(:,2)=x(:,2)*ls/l;
semilogy(x(:,1),x(:,2),'LineWidth',2,'DisplayName','Eps Eri (K2V, active)');

x=importdata('gj667c.txt');
l=trapz(x(:,1),x(:,2));
x(:,2)=x(:,2)*ls/l;
semilogy(x(:,1),x(:,2),'LineWidth',2,'DisplayName','GJ 667C (M1.5V)');

x=importdata('gj436.txt');
l=trapz(x(:,1),x(:,2));
x(:,2)=x(:,2)*ls/l;
semilogy(x(:,1),x(:,2),'LineWidth',2,'DisplayName','GJ 436 (M3.5V)');

x=importdata('ad_leo.txt');
l=trapz(x(:,1),x(:,2));
x(:,2)=x(:,2)*ls/l;
semilogy(x(:,1),x(:,2),'LineWidth',2,'DisplayName','AD Leo (M4.5V, active)');

x=importdata('gj876.txt');
l=trapz(x(:,1),x(:,2));
x(:,2)=x(:,2)*ls/l;
semilogy(x(:,1),x(:,2),'LineWidth',2,'DisplayName','GJ 876 (M5V)');



legend1 = legend(axes1,'show');
set(legend1,'EdgeColor',[1 1 1]);