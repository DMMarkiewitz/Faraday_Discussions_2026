close all
clear all

%%%%%%%%%%%%%%%%
%Key File Inputs
%%%%%%%%%%%%%%%%
mm = 21; %Molality
act_mesh = 'S_WiSE_Case_fp4_sticky_alpha_03_21m_01.mat'; %Mesh name
qsr = -0.2; % Surface Charge Density C/m^2
save_name = 'Sticky_WiSE_21m_01_a_03_m_02.mat'; %Activity Mesh save name
fig_save_name = 'Sticky_WiSE_21m_01_a_03_qs_m_02_act_all_no_ref';%Figure Save name

%%%%%%%%%%%%%%%%%%%
%Physical Constants
%%%%%%%%%%%%%%%%%%%
T = 300;         %temperature
Na = 6.022e23;   %Avogadro's number
kB = 1.38e-23;   %Boltzmann constant
e = 1.6e-19;     %elemntary charge
beta = 1/(kB*T); %inverse thermal energy

%General Systems Properties and Bulk Quantities
fa = 3;         %anion functionality
fc = 4;         %cation functionality
xic = 0.4;      %volume ratio of cation, relative to water, taken ot be 1 lattice 
xia = 10.8;     %volume ratio of anion, relative to water, taken ot be 1 lattice
MMw=18.02/1000; %kg/mol molar mass of water

%Property Selector
if mm == 12
    phi_p_b = 0.02526;
    v0 = 2.2863e-29; %m^3 
    L = 0.2259; %Lpm/Lpo; %Ion-Association Constant
elseif mm == 15
    phi_p_b = 0.02684;
    v0 = 2.2456e-29; %m^3 
    L = 0.2312; %Lpm/Lpo; %Ion-Association Constant
elseif mm == 21
    phi_p_b = 0.02888;
    v0 = 2.1744e-29; %m^3 
    L = 0.2527; %Lpm/Lpo; %Ion-Association Constant
else
    phi_p_b = 1/(1+(1/(mm*MMw) + xia)/xic);
    v0 = (2.1744e-29+2.2863e-29+2.2456e-29)/3; %m^3 volume
    L = (0.2259+0.2312+0.2527)/3; %Lpm/Lpo; %Ion-Association Constant
end

phi_m_b = xia/xic*phi_p_b;
psi_p_b = fc*phi_p_b/xic;
psi_m_b = fa*phi_m_b/xia;
phi_0_b = 1 - phi_p_b - phi_m_b;
phi_bulk = [phi_p_b,phi_m_b,phi_0_b];

%%%%%%%
%Inputs
%%%%%%%
eps_r = 10.1; %relative dielectric constant
eps=eps_r*8.85e-12; %full dielectric constant
%qsr = -0.2;
p = 1.67e-29/(sqrt(v0*eps/(beta*(phi_p_b/xic + phi_m_b/xia)))); %Dimensionless effected by varying the debye length, the dipole moment and the fundamental charge    
a=0.3; %charge rescaling parameter from Goodwin-Kornyshev model

%Mapping functions set-up
I_data = load(act_mesh);
[mPhi,mDPhi] = ndgrid(I_data.Phi,I_data.DPhi);
f_I_rho= griddedInterpolant(flip(mPhi),mDPhi,flip(I_data.rho),'makima');
f_I_phi_c = griddedInterpolant(flip(mPhi),mDPhi,flip(I_data.phi(:,:,1)),'makima');
f_I_phi_a = griddedInterpolant(flip(mPhi),mDPhi,flip(I_data.phi(:,:,2)),'makima'); 
f_I_phi_0 = griddedInterpolant(flip(mPhi),mDPhi,flip(I_data.phi(:,:,3)),'makima'); 
f_I_pc = griddedInterpolant(flip(mPhi),mDPhi,flip(I_data.pc_a),'makima');
f_I_pa = griddedInterpolant(flip(mPhi),mDPhi,flip(I_data.pa_a),'makima');
f_I_pp0 = griddedInterpolant(flip(mPhi),mDPhi,flip(I_data.eta_a),'makima'); 
f_I_p0 = griddedInterpolant(flip(mPhi),mDPhi,flip(I_data.p0_a),'makima'); 
f_I_epsr = griddedInterpolant(flip(mPhi),mDPhi,flip(I_data.epsr_a),'makima');
f_I_tau = griddedInterpolant(flip(mPhi),mDPhi,flip(I_data.tau_a),'makima');

%Accounting for symmetric nature of these function
f_I_rho = @(xPhi,xDPhi) f_I_rho(xPhi,abs(xDPhi));
f_I_phi_c = @(xPhi,xDPhi) f_I_phi_c(xPhi,abs(xDPhi));
f_I_phi_a = @(xPhi,xDPhi) f_I_phi_a(xPhi,abs(xDPhi));
f_I_phi_0 = @(xPhi,xDPhi) f_I_phi_0(xPhi,abs(xDPhi));
f_I_pc = @(xPhi,xDPhi) f_I_pc(xPhi,abs(xDPhi));
f_I_pa = @(xPhi,xDPhi) f_I_pa(xPhi,abs(xDPhi));
f_I_pp0 = @(xPhi,xDPhi) f_I_pp0(xPhi,abs(xDPhi));
f_I_p0 = @(xPhi,xDPhi) f_I_p0(xPhi,abs(xDPhi));
f_I_epsr = @(xPhi,xDPhi) f_I_epsr(xPhi,abs(xDPhi));
f_I_tau = @(xPhi,xDPhi) f_I_tau(xPhi,abs(xDPhi));

%Need to form secondary meshes
I_epsr = flip(I_data.epsr_a);
I_depsrdu = nan(size(I_data.eta_a));
I_depsrddu = nan(size(I_data.eta_a));
%Setting up order 6 stencils
fd_6 = [-49/20, 6, -15/2, 20/3, -15/4, 6/5, -1/6];
cen_6 = [-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60];
bd_6 = [1/6, -6/5, 15/4, -20/3, 15/2, -6, 49/20];
h_u = I_data.Phi(1)-I_data.Phi(2);
h_du = I_data.DPhi(2)-I_data.DPhi(1);

%Getting depsr/du mesh accuracy order 6
for idu = 1:length(I_data.DPhi)
    %LHS
    I_depsrdu(1:3,idu) = sum(fd_6.*[I_epsr(1:3,idu),I_epsr(2:4,idu),I_epsr(3:5,idu),I_epsr(4:6,idu),I_epsr(5:7,idu),I_epsr(6:8,idu),I_epsr(7:9,idu)],2)/h_u;
    %CEN
    I_depsrdu(4:end-3,idu) = sum(cen_6.*[I_epsr(1:end-6,idu),I_epsr(2:end-5,idu),I_epsr(3:end-4,idu),I_epsr(4:end-3,idu),I_epsr(5:end-2,idu),I_epsr(6:end-1,idu),I_epsr(7:end,idu)],2)/h_u;
    %RHS
    I_depsrdu(end-2:end,idu) = sum(bd_6.*[I_epsr(end-8:end-6,idu),I_epsr(end-7:end-5,idu),I_epsr(end-6:end-4,idu),I_epsr(end-5:end-3,idu),I_epsr(end-4:end-2,idu),I_epsr(end-3:end-1,idu),I_epsr(end-2:end,idu)],2)/h_u;
end

%Getting depsr/ddu mesh accuracy order 6
for iu = 1:length(I_data.Phi)
    %LHS
    I_depsrddu(iu,1:3) = sum(fd_6.*[I_epsr(iu,1:3);I_epsr(iu,2:4);I_epsr(iu,3:5);I_epsr(iu,4:6);I_epsr(iu,5:7);I_epsr(iu,6:8);I_epsr(iu,7:9)]',2)/h_du;
    %CEN
    I_depsrddu(iu,4:end-3) = sum(cen_6.*[I_epsr(iu,1:end-6);I_epsr(iu,2:end-5);I_epsr(iu,3:end-4);I_epsr(iu,4:end-3);I_epsr(iu,5:end-2);I_epsr(iu,6:end-1);I_epsr(iu,7:end)]',2)/h_du;
    %RHS
    I_depsrddu(iu,end-2:end) = sum(bd_6.*[I_epsr(iu,end-8:end-6);I_epsr(iu,end-7:end-5);I_epsr(iu,end-6:end-4);I_epsr(iu,end-5:end-3);I_epsr(iu,end-4:end-2);I_epsr(iu,end-3:end-1);I_epsr(iu,end-2:end)]',2)/h_du;
end

f_I_depsrdu = griddedInterpolant(flip(mPhi),mDPhi,I_depsrdu,'makima'); 
f_I_depsrddu = griddedInterpolant(flip(mPhi),mDPhi,I_depsrddu,'makima'); 
%Accounting for symmetric nature of these function
f_I_depsrdu = @(xPhi,xDPhi) f_I_depsrdu(xPhi,abs(xDPhi));
f_I_depsrddu = @(xPhi,xDPhi) sign(xDPhi)*f_I_depsrddu(xPhi,abs(xDPhi));


%Simulation Distance till we say we are in bulk
cutoff = ceil(2.1/sqrt(v0*eps/(beta*e^2*(phi_p_b/xic + phi_m_b/xia)))*10^-9);

%Step and Function tolerance should be same as grid solver
Tol = 1e-12;

%Runner Section
xGuess=[];
yGuess=[];

%Functions we need 
select_c = @(x,y,z) x{y}(z);

options = optimoptions('fsolve','display','off','TolX',Tol,'TolFun',Tol,'MaxFunctionEvaluations',1e6,'MaxIter',1e6);

%Solving for the qs equivalent to the qsr of interest
out_cell = @(x) cellify_MPB(@(x) MPB(xGuess,yGuess,cutoff,x,p,fc,fa,xic,xia,phi_bulk,f_I_rho,f_I_phi_c,f_I_phi_a,f_I_phi_0,f_I_pc,f_I_pa,f_I_pp0,f_I_p0,f_I_epsr,f_I_depsrdu,f_I_depsrddu,options),x);  %Building key element
if qsr>0
    [qs, ~, flag_qs, ~] = fzero(@(x) qsr-select_c(out_cell(x),7,1)*eps/(phi_p_b/xic + phi_m_b/xia)*x/sqrt(v0*eps/(beta*e^2*(phi_p_b/xic + phi_m_b/xia)))/beta/e,qsr);
elseif qsr<0
    [qs, ~, flag_qs, ~] = fzero(@(x) qsr-select_c(out_cell(x),7,1)*eps/(phi_p_b/xic + phi_m_b/xia)*x/sqrt(v0*eps/(beta*e^2*(phi_p_b/xic + phi_m_b/xia)))/beta/e,qsr);
else
    qs = 0;
    flag_qs = 1;
end

if flag_qs ~= 1
    warning('Not at true qs of interest')
    qs
end

[xint, Phi, DPhi, phi, prob, ppmpmp, epsr, rho] = MPB(xGuess,yGuess,cutoff,qs,p,fc,fa,xic,xia,phi_bulk,f_I_rho,f_I_phi_c,f_I_phi_a,f_I_phi_0,f_I_pc,f_I_pa,f_I_pp0,f_I_p0,f_I_epsr,f_I_depsrdu,f_I_depsrddu,options);

%We have to calculate the dimensional qs at the end
t = -epsr(1)*eps/(phi_p_b/xic + phi_m_b/xia)*DPhi(1)/sqrt(v0*eps/(beta*e^2*(phi_p_b/xic + phi_m_b/xia)))/beta/e;

dl = sqrt(v0*eps/(beta*e^2*(phi_p_b/xic + phi_m_b/xia)));% in m

%Building Saved Data
data.mu_p = 1+log((1+fc/xic)*phi(1,:).*(prob(3,:)).^fc)-fc*log(phi(3,:).*(1-prob(4,:)))-xic*log(f_I_tau(Phi,DPhi));
data.mu_m = 1+log(phi(2,:).*(1-prob(2,:)).^fa)-xia*log(f_I_tau(Phi,DPhi));
data.mu_0 = 1+log(phi(3,:).*(1-prob(4,:)))-log(f_I_tau(Phi,DPhi));
data.spatial_nm = xint*dl*10^9;%in nm

%Saving Activities for Merge Plotting
save(save_name,"data")

%Activity Predictions without reference
figure()
plot(xint*dl*10^9,1+log((1+fc/xic)*phi(1,:).*(prob(3,:)).^fc)-fc*log(phi(3,:).*(1-prob(4,:)))-xic*log(f_I_tau(Phi,DPhi)),'r-','linewidth',1.5)
hold on
plot(xint*dl*10^9,1+log(phi(2,:).*(1-prob(2,:)).^fa)-xia*log(f_I_tau(Phi,DPhi)),'b--','linewidth',1.5)
plot(xint*dl*10^9,1+log(phi(3,:).*(1-prob(4,:)))-log(f_I_tau(Phi,DPhi)),'k:','linewidth',1.5)
xlim([0,2])
temp = ylim;
ylim([temp(1)-.15,temp(2)+2.5])
box on
leg=legend({'Li$^+$','TFSI$^-$','H$_2$O'},'location','north','interpreter','latex','Orientation','horizontal','box','off');
leg.ItemTokenSize = [18,18];
pbaspect([(1+sqrt(5))/2,1,1])
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'DefaultTextInterpreter', 'latex');
xlabel(['Distance from electrode, nm'],'interpreter','latex', 'fontsize', 16)
ylabel(['Ln($\bar{a}_i$/$a_i^\theta$)'],'interpreter','latex', 'fontsize', 16)
g=gcf;
g.Renderer='painters';
set(gca, 'Xcolor', 'k');
set(gca, 'Ycolor', 'k');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1.5);
set(gca, 'Layer', 'Top');
set(gca, 'color', 'none');
set(g,'InvertHardcopy','on');
exportgraphics(gcf,[fig_save_name,'.eps'],'Resolution',600)
exportgraphics(gcf,[fig_save_name,'.jpeg'],'Resolution',600)
saveas(gcf,[fig_save_name,'.fig'])

function [xint, Phi, DPhi, phi, prob, ppmpmp, epsr, rho] = MPB(xGuess,yGuess,cutoff,qs,p,fc,fa,xic,xia,phi_bulk,f_I_rho,f_I_phi_c,f_I_phi_a,f_I_phi_0,f_I_pc,f_I_pa,f_I_pp0,f_I_p0,f_I_epsr,f_I_depsrdu,f_I_depsrddu,options)
    ode = @(x,u) odeMPB(x,u,p,fc,fa,xic,xia,phi_bulk,f_I_rho,f_I_phi_c,f_I_phi_a,f_I_phi_0,f_I_pc,f_I_pa,f_I_pp0,f_I_p0,f_I_epsr,f_I_depsrdu,f_I_depsrddu,options);
    bc =  @(ya,yb) bcMPB(ya,yb,qs,fc,fa,xic,xia,phi_bulk,f_I_rho,f_I_phi_c,f_I_phi_a,f_I_phi_0,f_I_pc,f_I_pa,f_I_pp0,f_I_p0,f_I_epsr,f_I_depsrdu,f_I_depsrddu,options);

    option=[];
    resol = cutoff*5;
    solinit = bvpinit(linspace(0,cutoff,resol),@initMPB, option, xGuess, yGuess, qs);
    sol = bvp4c(ode, bc, solinit, option);

    Inter_factor = 5;
    xint = linspace(0,cutoff,Inter_factor*resol);
    Sxint = deval(sol,xint);
    Phi = Sxint(1,:);
    DPhi = Sxint(2,:);
    epsr = zeros(1,length(Sxint(2,:)));
    for i = 1:length(Sxint(2,:))
        %So we can back calc qs with units of C/m^2
        epsr(i) = f_I_epsr(Phi(i),DPhi(i));
    end
    phi = nan(3,length(Phi)); % stores +,-,0 in different rows
    prob = nan(4,length(Phi)); % stores p+-,p-+,p+0,p0+ in different rows
    rho = nan(1,length(Sxint(1,:)));
    ppmpmp = nan(1,length(Sxint(1,:)));
    for i = 1:Inter_factor*resol
        P_out = [f_I_pc(Phi(i),DPhi(i)),f_I_pa(Phi(i),DPhi(i)),f_I_pp0(Phi(i),DPhi(i)),f_I_p0(Phi(i),DPhi(i))];
        Phi_out = [f_I_phi_c(Phi(i),DPhi(i)),f_I_phi_a(Phi(i),DPhi(i)),f_I_phi_0(Phi(i),DPhi(i))];

        phi(1,i) = Phi_out(1);
        phi(2,i) = Phi_out(2);
        phi(3,i) = Phi_out(3);

        prob(1,i) = 1-P_out(3);
        prob(2,i) = P_out(2);
        prob(3,i) = P_out(3);
        prob(4,i) = P_out(4);

        rho(i) = f_I_rho(Phi(i),DPhi(i));

        ppmpmp(i) =  prob(1,i)*prob(2,i);
        
    end
    
end

function dudx = odeMPB(x,u,p,fc,fa,xic,xia,phi_bulk,f_I_rho,f_I_phi_c,f_I_phi_a,f_I_phi_0,f_I_pc,f_I_pa,f_I_pp0,f_I_p0,f_I_epsr,f_I_depsrdu,f_I_depsrddu,options)
    P_out = [f_I_pc(u(1),u(2)),f_I_pa(u(1),u(2)),f_I_pp0(u(1),u(2)),f_I_p0(u(1),u(2))];
    Phi_out = [f_I_phi_c(u(1),u(2)),f_I_phi_a(u(1),u(2)),f_I_phi_0(u(1),u(2))];

    phi_c = Phi_out(1);
    phi_a = Phi_out(2);
    phi_0 = Phi_out(3); 

    pce = 1-P_out(3);
    pae = P_out(2);
        
    dudx=[u(2)
          -((f_I_rho(u(1),u(2))) + u(2).^2 * f_I_depsrdu(u(1),u(2)))/(f_I_depsrddu(u(1),u(2))*u(2) + f_I_epsr(u(1),u(2)))];
end

function res = bcMPB(ya,yb,qs,fc,fa,xic,xia,phi_bulk,f_I_rho,f_I_phi_c,f_I_phi_a,f_I_phi_0,f_I_pc,f_I_pa,f_I_pp0,f_I_p0,f_I_epsr,f_I_depsrdu,f_I_depsrddu,options)
    if qs~=0
        res = [ya(2)+qs %charge of surface
               yb(1)]; %far away
    else
        res = [ya(2) %charge of surface
        yb(1)]; %far away
    end

end

function yinit = initMPB(x, xGuess, yGuess,qs)
    if length(xGuess)>0
        yinit = pchip(xGuess, yGuess,x);
    else
        yinit = [qs*exp(-x)
                 -qs*exp(-x)];
    end
end

function out_cell = cellify_MPB(fun,x)
    %Packs the function handle of MPB as a function of qs into cells
    [xint, Phi, DPhi, phi, prob, ppmpmp, epsr, rho] = fun(x);
    out_cell = {xint, Phi, DPhi, phi, prob, ppmpmp, epsr, rho};
end