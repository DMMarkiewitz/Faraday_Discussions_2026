function Implicit_Poisson_Solver
    clear all
    close all

    %Key File Inputs
    mm = 21; %Molality
    act_mesh = 'S_WiSE_Case_fp4_sticky_alpha_03_21m_01.mat'; %Mesh name
    qsr = 0.2; % Surface Charge Density C/m^2
    fig_save_name = 'New_Final_WiSE_21m_fp4_sticky_cl_dist_p_';%Figure Save name
    fig_save_name_qsr = '2'; %10*qsr
    fig_save_name_mesh_size = '01'; %grid spacing
    
    %Physical Constants
    T = 300; 
    Na = 6.022e23; 
    kB = 1.38e-23;  
    e = 1.6e-19;
    beta = 1/(kB*T);
    EP = 50;
    Ep_flag = 0;

    %General Systems Properties and Bulk Quantities
    fa = 3;
    fc = 4;
    xic = 0.4;
    xia = 10.8;
    
    if mm == 12
        phi_p_b = 0.02526;
        v0 = 2.2863e-29; %m^3 volume backcalculated from vbox and Omega
        L = 0.2259; %Lpm/Lpo; %Ion-Association Constant
    elseif mm == 15
        phi_p_b = 0.02684;
        v0 = 2.2456e-29; %m^3 volume backcalculated from vbox and Omega
        L = 0.2312; %Lpm/Lpo; %Ion-Association Constant
    elseif mm == 21
        phi_p_b = 0.02888;
        v0 = 2.1744e-29; %m^3 volume backcalculated from vbox and Omega
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
    
    %Inputs
    eps_r = 10.1;
    eps=eps_r*8.85e-12;
    %qsr = 0.2;
    p = 1.67e-29/(sqrt(v0*eps/(beta*(phi_p_b/xic + phi_m_b/xia)))); %Dimensionless effected by varying the debye length, the dipole moment and the fundamental charge    
    a=0.3;
    spec_pt = 5/sqrt(v0*eps/(beta*e^2*(phi_p_b/xic + phi_m_b/xia)))*10^-10;

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
    cutoff = 30; 

    %Step and Function tolerance should be same as grid solver
    Tol = 1e-12;
    
    %Runner Section%
    xGuess=[];
    yGuess=[];

    %Functions we need 
    select_c = @(x,y,z) x{y}(z);
   
    options = optimoptions('fsolve','display','off','TolX',Tol,'TolFun',Tol,'MaxFunctionEvaluations',1e6,'MaxIter',1e6);

    %Solving for the qs equivalent to the qsr of interest
    out_cell = @(x) cellify_MPB(@(x) MPB(xGuess,yGuess,cutoff,x,p,L,fc,fa,xic,xia,phi_bulk,f_I_rho,f_I_phi_c,f_I_phi_a,f_I_phi_0,f_I_pc,f_I_pa,f_I_pp0,f_I_p0,f_I_epsr,f_I_depsrdu,f_I_depsrddu,v0,Ep_flag,EP,spec_pt,options),x);  %Building key element
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

    %Turning on EP flag
    Ep_flag = 0;
    [xint, Phi, DPhi, phi, prob, ppmpmp, epsr, rho, AggLs_i, clms_a_i] = MPB(xGuess,yGuess,cutoff,qs,p,L,fc,fa,xic,xia,phi_bulk,f_I_rho,f_I_phi_c,f_I_phi_a,f_I_phi_0,f_I_pc,f_I_pa,f_I_pp0,f_I_p0,f_I_epsr,f_I_depsrdu,f_I_depsrddu,v0,Ep_flag,EP/10,spec_pt,options);

    %We have to calculate the dimensional qs at the end
    t = -epsr(1)*eps/(phi_p_b/xic + phi_m_b/xia)*DPhi(1)/sqrt(v0*eps/(beta*e^2*(phi_p_b/xic + phi_m_b/xia)))/beta/e

    dl = sqrt(v0*eps/(beta*e^2*(phi_p_b/xic + phi_m_b/xia)))*10^10;% in A

    %Three points for cutting
    [~,cI1] = min(abs(2-xint));
    [~,cI2] = min(abs(5-xint));
    [~,cI3] = min(abs(10-xint));%7.5
    [~,cInt] = min(abs(5-xint*dl));%7.5
    CI = [cI1,cI2,cI3];
    for i_c = 1:3 
        Lbar = L * sinch(p*DPhi(CI(i_c)));
        zp = phi(3,CI(i_c)) * (1-prob(4,CI(i_c))) / Lbar;
        Xp = fc/xic * phi(1,CI(i_c)) * (1-prob(1,CI(i_c)))^fc / zp;
        Yp = fa/xia * phi(2,CI(i_c)) * (1-prob(2,CI(i_c)))^fa / zp;
        [clms,cubeAggLsdv0] = WiSE_sti_theory_func_count(fa,fc,xic,xia,Lbar,Xp,Yp,zp,EP);
        AggLs(i_c) = (cubeAggLsdv0*v0)^(1/3);
        clms_a{i_c} = clms;
    end
    %Obtaining the mean over 5A
    clms_5 = zeros(51,51);
    for i=1:cInt
        Lbar = L * sinch(p*DPhi(i));
        zp = phi(3,i) * (1-prob(4,i)) / Lbar;
        Xp = fc/xic * phi(1,i) * (1-prob(1,i))^fc / zp;
        Yp = fa/xia * phi(2,i) * (1-prob(2,i))^fa / zp;
        [clms,cubeAggLsdv0] = WiSE_sti_theory_func_count(fa,fc,xic,xia,Lbar,Xp,Yp,zp,EP);
        clms_5 = clms_5 + clms/length(1:cInt);%clms_a_i{i};
    end
    %clms_5 = mean(clms_i5,3);

    %Plotting averaged CDF
    figure()
    temp=log10(clms_5);%log10(clms_a{cI});
    temp(temp<-8)=nan;%8%Need to choose our cutoff magnitude
    fig=pcolor(-0.5:size(temp,1)-1.5,-0.5:size(temp,2)-1.5,temp')
    set(fig, 'EdgeColor', 'none')
    clim([-8,0])
    xlim([0,20])
    ylim([0,20])
    cb = colorbar("LineWidth",1.5)
    ylabel(cb,'Log$_{10}$($c_{lms}$)', 'interpreter','latex', 'fontsize', 16)
    hold on
    plot([0,100],[0,100],'--k','LineWidth',2)
    ylabel('number of anions in cluster, $m$', 'interpreter','latex', 'fontsize', 16)
    xlabel('number of cations in cluster, $l$', 'interpreter','latex', 'fontsize', 16)
    %Can add more
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    box on
    ax = gca;
    g=gcf;
    g.Renderer='painters';
    ax.FontSize = 16;
    %title('Theory', 'interpreter','latex', 'fontsize', titlesize)
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'color', 'none');
    set(g,'InvertHardcopy','on');

    exportgraphics(g,[fig_save_name,fig_save_name_qsr,'_avg_5_',fig_save_name_mesh_size,'.eps'],'Resolution',600)
    exportgraphics(g,[fig_save_name,fig_save_name_qsr,'_avg_5_',fig_save_name_mesh_size,'.jpeg'],'Resolution',600)
    saveas(g,[fig_save_name,fig_save_name_qsr,'_avg_5_',fig_save_name_mesh_size,'.fig'])

    % Importing / Defining Experimental Data
    %[~,cI] = min(abs(2-xint))
    figure()
    temp=log10(clms_a{1})%log10(clms_a{cI});
    temp(temp<-8)=nan;%8%Need to choose our cutoff magnitude
    fig=pcolor(-0.5:size(temp,1)-1.5,-0.5:size(temp,2)-1.5,temp')
    set(fig, 'EdgeColor', 'none')
    clim([-8,0])
    xlim([0,20])
    ylim([0,20])
    cb = colorbar("LineWidth",1.5)
    ylabel(cb,'Log$_{10}$($c_{lms}$)', 'interpreter','latex', 'fontsize', 16)
    hold on
    plot([0,100],[0,100],'--k','LineWidth',2)
    ylabel('number of anions in cluster, $m$', 'interpreter','latex', 'fontsize', 16)
    xlabel('number of cations in cluster, $l$', 'interpreter','latex', 'fontsize', 16)
    %Can add more
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    box on
    ax = gca;
    g=gcf;
    g.Renderer='painters';
    ax.FontSize = 16;
    %title('Theory', 'interpreter','latex', 'fontsize', titlesize)
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'color', 'none');
    set(g,'InvertHardcopy','on');

    exportgraphics(g,[fig_save_name,fig_save_name_qsr,'_place_2_',fig_save_name_mesh_size,'.eps'],'Resolution',600)
    exportgraphics(g,[fig_save_name,fig_save_name_qsr,'_place_2_',fig_save_name_mesh_size,'.jpeg'],'Resolution',600)
    saveas(g,[fig_save_name,fig_save_name_qsr,'_place_2_',fig_save_name_mesh_size,'.fig'])

    % Importing / Defining Experimental Data
    %[~,cI] = min(abs(5-xint))
    figure()
    temp=log10(clms_a{2});%log10(clms_a{cI});
    temp(temp<-8)=nan;%8%Need to choose our cutoff magnitude
    fig=pcolor(-0.5:size(temp,1)-1.5,-0.5:size(temp,2)-1.5,temp')
    set(fig, 'EdgeColor', 'none')
    clim([-8,0])
    xlim([0,20])
    ylim([0,20])
    cb = colorbar("LineWidth",1.5)
    ylabel(cb,'Log$_{10}$($c_{lms}$)', 'interpreter','latex', 'fontsize', 16)
    hold on
    plot([0,100],[0,100],'--k','LineWidth',2)
    ylabel('number of anions in cluster, $m$', 'interpreter','latex', 'fontsize', 16)
    xlabel('number of cations in cluster, $l$', 'interpreter','latex', 'fontsize', 16)
    %Can add more
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    box on
    ax = gca;
    g=gcf;
    g.Renderer='painters';
    ax.FontSize = 16;
    %title('Theory', 'interpreter','latex', 'fontsize', titlesize)
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'color', 'none');
    set(g,'InvertHardcopy','on');

    exportgraphics(g,[fig_save_name,fig_save_name_qsr,'_place_5_',fig_save_name_mesh_size,'.eps'],'Resolution',600)
    exportgraphics(g,[fig_save_name,fig_save_name_qsr,'_place_5_',fig_save_name_mesh_size,'.jpeg'],'Resolution',600)
    saveas(g,[fig_save_name,fig_save_name_qsr,'_place_5_',fig_save_name_mesh_size,'.fig'])

    % Importing / Defining Experimental Data
    figure()
    temp=log10(clms_a{3});
    temp(temp<-8)=nan;
    fig=pcolor(-0.5:size(temp,1)-1.5,-0.5:size(temp,2)-1.5,temp')
    set(fig, 'EdgeColor', 'none')
    clim([-8,0])
    xlim([0,20])
    ylim([0,20])
    cb = colorbar("LineWidth",1.5)
    ylabel(cb,'Log$_{10}$($c_{lms}$)', 'interpreter','latex', 'fontsize', 16)
    hold on
    plot([0,100],[0,100],'--k','LineWidth',2)
    ylabel('number of anions in cluster, $m$', 'interpreter','latex', 'fontsize', 16)
    xlabel('number of cations in cluster, $l$', 'interpreter','latex', 'fontsize', 16)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    box on
    ax = gca;
    g=gcf;
    g.Renderer='painters';
    ax.FontSize = 16;
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'color', 'none');
    set(g,'InvertHardcopy','on');

    exportgraphics(g,[fig_save_name,fig_save_name_qsr,'_place_10_',fig_save_name_mesh_size,'.eps'],'Resolution',600)
    exportgraphics(g,[fig_save_name,fig_save_name_qsr,'_place_10_',fig_save_name_mesh_size,'.jpeg'],'Resolution',600)
    saveas(g,[fig_save_name,fig_save_name_qsr,'_place_10_',fig_save_name_mesh_size,'.fig'])

end

function [xint, Phi, DPhi, phi, prob, ppmpmp, epsr, rho, AggLs, clms_a] = MPB(xGuess,yGuess,cutoff,qs,p,L,fc,fa,xic,xia,phi_bulk,f_I_rho,f_I_phi_c,f_I_phi_a,f_I_phi_0,f_I_pc,f_I_pa,f_I_pp0,f_I_p0,f_I_epsr,f_I_depsrdu,f_I_depsrddu,v0,Ep_flag,EP,spec_pt,options)
    ode = @(x,u) odeMPB(x,u,p,fc,fa,xic,xia,phi_bulk,f_I_rho,f_I_phi_c,f_I_phi_a,f_I_phi_0,f_I_pc,f_I_pa,f_I_pp0,f_I_p0,f_I_epsr,f_I_depsrdu,f_I_depsrddu,options);
    bc =  @(ya,yb) bcMPB(ya,yb,qs,fc,fa,xic,xia,phi_bulk,f_I_rho,f_I_phi_c,f_I_phi_a,f_I_phi_0,f_I_pc,f_I_pa,f_I_pp0,f_I_p0,f_I_epsr,f_I_depsrdu,f_I_depsrddu,options);

    option=[];
    resol = cutoff*5;
    solinit = bvpinit(linspace(0,cutoff,resol),@initMPB, option, xGuess, yGuess, qs);
    sol = bvp4c(ode, bc, solinit, option);

    Inter_factor = 5;
    xint = linspace(0,cutoff,Inter_factor*resol);
    xint = sort([xint,spec_pt]);
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
    AggLs = nan(1,length(Sxint(1,:)));
    clms_a = cell(1,length(Sxint(1,:)));

    for i = 1:Inter_factor*resol+1
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

        if Ep_flag == 1
            Lbar = L * sinch(p*DPhi(i));
            zp = phi(3,i) * (1-prob(4,i)) / Lbar;
            Xp = fc/xic * phi(1,i) * (1-prob(1,i))^fc / zp;
            Yp = fa/xia * phi(2,i) * (1-prob(2,i))^fa / zp;
            [clms,cubeAggLsdv0] = WiSE_sti_theory_func_count(fa,fc,xic,xia,Lbar,Xp,Yp,zp,EP);
            AggLs(i) = (cubeAggLsdv0*v0)^(1/3);
            clms_a{i} = clms;
        end
        
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

function B = sinch(x)
    B=x;
    for i=1:length(x)
        if abs(x(i)) > 1e-3
            B(i) = sinh(x(i))./x(i);
        else
            B(i) = 1 + x(i).^2 / 6 + x(i).^4 / 120 + x(i).^6 / 5040;
        end
    end
end

function [clms,cubeAggLs] = WiSE_sti_theory_func_count(fa,fc,xic,xia,Lbar,Xp,Yp,zp,EP)
clms = zeros(EP+1,EP+1);
cubeAggLs = 0;
cl = zeros([EP EP]);
cm = zeros([EP EP]);
cs = zeros([EP EP]);
    for J = 0:EP
        for K = 0:EP
            l = J+K;
            if l > 0 && fc * J - J - K + 1 >= 0 && fa * K - K - J + 1 >=0
                Temp = zp*sti_stockmeyer_prob_4sum(J,K,Xp,Yp,fc,fa);

                clms(J+1,K+1) = Temp;
                cubeAggLs = cubeAggLs + (xic*J+xia*K+fc*J-J-K+1)^2 * Temp;
                cl(J+1,K+1) = J*Temp;
                cm(J+1,K+1) = K*Temp;
                cs(J+1,K+1) = (fc*J-J-K+1)*Temp;

            elseif J == 0 && K == 0
                cs(J+1,K+1) = zp*Lbar;
                clms(J+1,K+1) = zp*Lbar;
                cubeAggLs = cubeAggLs + zp*Lbar;
            end
        end
    end
end

function SN = sti_stockmeyer_prob_4sum(l,m,XX,YY,fc,fa)

    SN = factorial(fa * m - m) .* factorial(fc * l - l)...
        ./factorial(fa * m - m - l + 1)./factorial(fc * l - l - m + 1)...
        * XX .^ l / factorial(l) * YY .^ m / factorial(m);

    if isnan(SN) || isinf(SN)

       lSN = lFactorial(fa * m - m)+lFactorial(fc * l - l)...
           -lFactorial(fa * m - m - l + 1)-lFactorial(fc * l - l - m + 1)...
           +l*log(XX)+m*log(YY)-lFactorial(l)-lFactorial(m);
       SN = exp(lSN);

    end
    
end

function A = lFactorial(n)

    if n > 170 
        A = sum(log(1:n));
    else
        A = log(factorial(n));
    end
    
end

function out_cell = cellify_MPB(fun,x)
    %Packs the function handle of MPB as a function of qs into cells
    [xint, Phi, DPhi, phi, prob, ppmpmp, epsr, rho, AggLs] = fun(x);
    out_cell = {xint, Phi, DPhi, phi, prob, ppmpmp, epsr, rho, AggLs};
end