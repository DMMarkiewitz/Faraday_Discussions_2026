%%%%%%%%%%%%%%%%%%%%
%Defining parameters
%%%%%%%%%%%%%%%%%%%%
clear all 
close all
clc
T = 300; 
Na = 6.022e23; 
kB = 1.38e-23;  
e = 1.6e-19;
beta = 1/(kB*T);
fa = 3;%3
fc = 4;%4
xic = 0.4;
xia = 10.8;
a = 0.4; %the correction for short range repulsion - going to set it to 1 to cause issues
eps_r = 10.1;
eps=eps_r*8.85e-12;
MMw=18.02/1000 %kg/mol molar mass of water

%range of search values
xv = 12;
for m=xv
    if m == 12
        phi_p_b = 0.02526;
        v0 = 2.2863e-29; %m^3 volume backcalculated from vbox and Omega,
        L = 0.2259 %Lpm/Lpo;
    elseif m == 15
        phi_p_b = 0.02684;
        v0 = 2.2456e-29; %m^3 volume backcalculated from vbox and Omega,
        L = 0.2312 %Lpm/Lpo; 
    elseif m == 21
        phi_p_b = 0.02888;
        v0 = 2.1744e-29; %m^3 volume backcalculated from vbox and Omega,
        L = 0.2527  %Lpm/Lpo;
    else
        phi_p_b = 1/(1+(1/(m*MMw) + xia)/xic);
        v0 = (2.1744e-29+2.2863e-29+2.2456e-29)/3; %m^3 volume backcalculated from vbox and Omega,
        L = (0.2259+0.2312+0.2527)/3; %Lpm/Lpo;
    end

    %Seeting up key values for the loop
    phi_m_b = xia*phi_p_b/xic;
    psi_p_b = fc*phi_p_b/xic;
    psi_m_b = fa*phi_m_b/xia;
    phi_0_b = 1 - phi_p_b - phi_m_b;%2*phi_s_b;
    p = 1.67e-29/(sqrt(v0*eps/(beta*(phi_p_b/xic + phi_m_b/xia))));
    
    %need to have Phi in increasing values
    %Must be symmetric in the new homopoty case else need to modify more
    re=0.75*e*beta;
    xi=-0.01;%-.01;
    Phi = re:xi:-0.75*e*beta;% 2:-0.01:-2;%-0.01:-0.01:-3, 
    DPhi = 0:0.01:2.5;%:0.01:0.01;%for now, going to ignore this....
    %%0:0.01:.8%8;%3.5;%0:0.01:2;
    l=length(Phi);
    %Storage matrix for the volume fractions
    phi = nan(length(Phi),length(DPhi),3);
    TAU = nan(length(Phi),length(DPhi));
    %Phi = -5;
    %Is EP just total number of ions in cluster used?
    EP = 500;%500 is okay, so is 250 for pre-gel, but smaller values could be used - does not matter much here!
    %volumes are also the same
    
    %Storage matrix for the eta here fractions
    eta_a = nan(length(Phi),length(DPhi),1);
    %Storage matrix for tau fractions
    tau_a = nan(length(Phi),length(DPhi),1);
    %Storing epsilon data
    epsr_a = nan(length(Phi),length(DPhi),1);
    %Storing pa data for testing
    pa_a = nan(length(Phi),length(DPhi),1);
    %Storing p0c data for testing
    p0_a = nan(length(Phi),length(DPhi),1);
    %Saving rho
    rho = nan(length(Phi),length(DPhi),1);
    
    %Step and Function tolerance
    %Tol = 1e-10;
    Tol_v = 1e-5;%1e-10;
    %used to insure no points fail to solve with current method
    T_flag = 0;
    %used to check that the sticky ion approximation is still valid
    Con_flag = 0;
    %used to check if gelation has occured
    Gel_flag = 0;
    %Gel and fails to have incompressibility solution
    Gel_np_flag = 0;
    %We fail to have solution with incompressibility
    F_flag = 0;
    
    %reindexing properly to allow for quicker evaluations through homotopy
    %and not needing two different vectors and combining them and such
    h_Phi = 1:l;
    h_Phi(1:length(re:xi:xi)) = flip(h_Phi(1:length(re:xi:xi)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Setting up the bulk properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pc = f_pc(phi_0_b,psi_p_b,psi_m_b,L);
    pa = f_pa(phi_0_b,psi_p_b,psi_m_b,L);
    pc0 = f_pc0(phi_0_b,psi_p_b,psi_m_b,L);
    p0 = f_p0(phi_0_b,psi_p_b,psi_m_b,L);
    P_temp = [pc0,pa,p0];%Defining initial prob vect
    f_g = pc0; %Initial initial guess
    f_g_0 = pc0; %Initial initial guess backup, used for continuation on individual branches


    %Gelation criteria
    1 - pa*pc*(fc-1)*(fa-1)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Finding charge density as a function of field
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ig=[1,phi_p_b,phi_m_b,phi_0_b];%we could give it a much better initial guess ...
    Test = 0;
    Tau=1;

    
    select = @(x,y) x(y);

    optionsfs = optimoptions('fsolve','display','off','TolX',Tol_v,'TolFun',Tol_v,'MaxFunctionEvaluations',1e4,'MaxIter',1e4);


    for j=1:length(DPhi) 
        j
        
        Lbar = L * sinch(p*DPhi(j))
       
        %Defining useful equation solver for precision problems
        syms tL_b(x,y,z)
        tL_b(x,y,z) = z - x*(1-y)/(y*(1-x));
        for i=h_Phi %Using list for loop to achieve dual directional continuaion homotopy
            %i
            %Phi(i)
            
            res_t = @(key) [key(2) + key(3) + key(4) - 1;
                            key(2) - phi_p_b*(1-pc)^fc*exp(-a*Phi(i))*key(1)^(xic+fc)/f_pc0(key(4),fc*key(2)/xic,fa*key(3)/xia,Lbar)^fc;
                            %key(2) - phi_p_b*(phi_m_b*f_pc(key(4),fc*key(2)/xic,fa*key(3)/xia,Lbar)*(1-pc))^fc*exp(-a*Phi(i))*key(1)^(xic)/(key(3)*pc*(1-f_pc(key(4),fc*key(2)/xic,fa*key(3)/xia,Lbar)))^fc;%this is all I would need to change
                            key(3) - phi_m_b*(1-pa)^fa*exp(a*Phi(i))*key(1)^(xia)/(1-f_pa(key(4),fc*key(2)/xic,fa*key(3)/xia,Lbar))^fa;
                            key(4) - sinch(p*DPhi(j))*phi_0_b*(1 - p0)*key(1)/(1-f_p0(key(4),fc*key(2)/xic,fa*key(3)/xia,Lbar))];
            
            %change the initial guess here?
            %IG = ig;
            if i == 2900 && j == 1
                IG = ig;
            end

            if i < 2900 && j == 1
                IG = [TAU(i+1,j),phi(i+1,j,1),phi(i+1,j,2),phi(i+1,j,3)];
            end

            if i > 2900 && j == 1
                IG = [TAU(i-1,j),phi(i-1,j,1),phi(i-1,j,2),phi(i-1,j,3)];
            end
            
            
            if j > 1
                if i <= 2900
                    IG = [TAU(i+1,j-1),phi(i+1,j-1,1),phi(i+1,j-1,2),phi(i+1,j-1,3)];
                end
    
                if i > 2900
                    IG = [TAU(i-1,j-1),phi(i-1,j-1,1),phi(i-1,j-1,2),phi(i-1,j-1,3)];
                end
            end

            %IG
            [keys,~,flag,~] = fsolve(@(x) res_t(x),abs(real(IG)),optionsfs)
            %if flag <= 0 %Sometimes it needs a bump in the right direction 
            %    [keys,~,flag,~] = fsolve(@(x) res_t(x),IG/2,optionsfs);
            %end
            %keys
            if flag > 0

                P_temp = [f_pc0(keys(4),fc*keys(2)/xic,fa*keys(3)/xia,Lbar),f_pa(keys(4),fc*keys(2)/xic,fa*keys(3)/xia,Lbar),f_p0(keys(4),fc*keys(2)/xic,fa*keys(3)/xia,Lbar)];

                phi_c = abs(real(keys(2)));
                %Phi(i)
                phi_0 = abs(real(keys(4)));
                phi_a = abs(real(keys(3)));

                %Defining phis
                phi(i,j,1) = phi_c;    
                phi(i,j,2) = phi_a;
                phi(i,j,3) = phi_0;
                TAU(i,j) = abs(real(keys(1)));

                %calculating rho
                rho(i,j) = (phi_c/xic - phi_a/xia);

                %Getting Probability
                pce = 1-P_temp(1);
                pae = P_temp(2);
                pa_a(i,j) = pae;
                p0_a(i,j) = P_temp(3);

                %Storing all the values
                eta_a(i,j) = P_temp(1);
                tau_a(i,j) = abs(real(keys(1)));
                epsr_a(i,j) = f_t_epsr(phi_p_b/xic,phi_m_b/xia,Phi(i),DPhi(j),p,phi_0,P_temp(3),select);

                Tau = keys(1);
                %abs(pae*(1-P_temp(3))/(P_temp(3)*(1-pae)))
                %Testing all original equations hold
                 %abs(fa*phi_a*pae/xia-fc*phi_c*pce/xic)
                 %abs(fc*phi_c*P_temp(1)/xic-phi_0*P_temp(3))
                 %abs(phi_a-(Tau^(xia)*phi_m_b*(1-pa)^fa *exp(a*Phi(i))/(1-P_temp(2))^fa))
                 %abs(phi_0 - (Tau*(phi_0_b*(1 - p0)*sinch(p*DPhi(j)))/(1-P_temp(3))))
                 %abs(1-phi_c-phi_a-phi_0)
                if abs(fa*phi_a*pae/xia-fc*phi_c*pce/xic) > Tol_v ||  abs(fc*phi_c*P_temp(1)/xic-phi_0*P_temp(3)) > Tol_v || abs(double(tL_b(pae,P_temp(3),Lbar))) > Tol_v
                    abs(fa*phi_a*pae/xia-fc*phi_c*pce/xic)
                    abs(fc*phi_c*P_temp(1)/xic-phi_0*P_temp(3))
                    abs(double(tL_b(pae,P_temp(3),Lbar)))
                    
                    Test = Test + 1;
                end
                
                %how come the cation is not checked here? because of
                %incompressibility?
                if abs(phi_a-(Tau^(xia)*phi_m_b*(1-pa)^fa *exp(a*Phi(i))/(1-P_temp(2))^fa)) > Tol_v || abs(phi_0 - (Tau*(phi_0_b*(1 - p0)*sinch(p*DPhi(j)))/(1-P_temp(3)))) > Tol_v || abs(1-phi_c-phi_a-phi_0) > Tol_v
                    Test = Test + sqrt(-1);
                end

                %Test to see if any outputs are nonphysical (Test element)
                if isnan(P_temp)
                    T_flag = T_flag + 1;
                    phi(i,j,1) = nan;    
                    phi(i,j,2) = nan;
                    phi(i,j,3) = nan;
                    eta_a(i,j) = nan;
                    tau_a(i,j) = nan;
                    epsr_a(i,j) = nan;
                elseif fc*phi_c/xic - fa*phi_a/xia - phi_0 > 0
                    Con_flag = Con_flag + 1;
                    phi(i,j,1) = nan;    
                    phi(i,j,2) = nan;
                    phi(i,j,3) = nan;
                    eta_a(i,j) = nan;
                    tau_a(i,j) = nan;
                    epsr_a(i,j) = nan;
                elseif phi_c+Tol_v < 0 || phi_0+Tol_v < 0 || phi_a+Tol_v < 0 || phi_a-Tol_v > 1 || phi_0-Tol_v > 1 || phi_c-Tol_v > 1 || ~isreal(phi_c) || ~isreal(phi_0) || ~isreal(phi_a)
                    T_flag = T_flag + 1;
                    phi(i,j,1) = nan;    
                    phi(i,j,2) = nan;
                    phi(i,j,3) = nan;
                    eta_a(i,j) = nan;
                    tau_a(i,j) = nan;
                    epsr_a(i,j) = nan;
                elseif (fa-1)*(fc-1)*pae*pce-1 > 0
                    Gel_flag = Gel_flag + 1;
                else
                    %f_g = p_out(Tau);
                    %ig = Tau;
                    
                    %ig = keys;
                    
                    % zp = phi_0 * (1-p0_a(i,j)) / Lbar;
                    % Xp = fc/xic * phi_c * (1-pce)^fc / zp;
                    % Yp = fa/xia * phi_a * (1-pae)^fa / zp;
                    % [~,ip_vf,agg_vf,sum_cl,sum_cm,sum_cs] = WiSE_sti_theory_func(fa,fc,Lbar,Xp,Yp,zp,EP)
                    % t=1
                end

                %New added to fix rare bug, if our guess is negative can cause
                %issues so making H_key always pos.
                %ig = abs(ig);

                %if i == length(re:xi:xi) %So we can smoothly solve our system, grabbing ref point
%                         f_g = p_out(Tau);
%                         f_g_0 = p_out(Tau);
                %    ig = keys;
                %    ig_0 = keys;
%                         P_temp = f_p_out(ig);
                %elseif i == 1%using ref point to staart second branch
%                         f_g = f_g_0;
                %    ig = ig_0;
%                         P_temp = f_p_out(ig);
                %elseif i == length(h_Phi) %restoring ref point for next branch to start correctly
%                         f_g = f_g_0;
                %    ig = ig_0;
                    %P_temp = f_p_out(ig);
                %end

            else

                i
                j
                Phi(i)
                error('Outside of solvable region')
                
            end
            
        end
    end
    
    Test
    
    T_flag
    
    Con_flag
    
    Gel_flag
        
    %Getting +-
    pc_a = 1-eta_a;
    
    %Saving the required
    save('WiSE_Case_fp4_sticky_alpha_0.4_12m_01_DM_2.mat','Phi','DPhi','rho','phi','pc_a','eta_a','pa_a','p0_a','tau_a','epsr_a')

    %Plotting figures for our system
    %Contours are hard set
    f_plot(Phi,DPhi,rho,phi,eta_a,pa_a,p0_a,tau_a,epsr_a,1-(fc-1)*(fa-1)*((1-eta_a).*pa_a))

end

function f_plot(Phi,DPhi,rho,phi,eta_a,pa_a,p0_a,tau_a,epsr_a,prox_gel)
    [Phi,DPhi] = meshgrid(Phi,DPhi);

    figure
    fig=gcf;
%     fig.Position(3) = 400;
%     fig.Position(4) = 400;
    contourf(Phi,DPhi,rho',250,'edgecolor','none')
    %plot(Phi,phi(:,:,3)','linewidth',1.5)
    hold on
    %colorbar
    %title('$\bar{\phi}_0$', 'interpreter','latex', 'fontsize', 16)
    %title('$\phi_0$ dependence on $\tilde{\nabla}u$ \& $u$', 'interpreter','latex', 'fontsize', 16)
    ylabel('$\tilde{\nabla}u$', 'interpreter','latex', 'fontsize', 16)
    xlabel('$u$', 'interpreter','latex', 'fontsize', 16)
    cb=colorbar();
    ylabel(cb,'$\tilde{\rho}$', 'interpreter','latex', 'fontsize', 16)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultTextInterpreter', 'latex')
    ax = gca;
    ax.FontSize = 16;
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'Units','normalized', 'Position',[0.1 0.15 0.675 0.8]);
    set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','on');
    % saveas(gcf,'Sticky_WiSE_05m_rho_Mesh','jpeg')
    % saveas(gcf,'Sticky_WiSE_05m_rho_Mesh','epsc')
    % saveas(gcf,'Sticky_WiSE_05m_rho_Mesh','fig')
    exportgraphics(gcf,['Sticky_WiSE_05m_rho_Mesh_01_a_01.eps'],'Resolution',600)
    exportgraphics(gcf,['Sticky_WiSE_05m_rho_Mesh_01_a_01.jpeg'],'Resolution',600)
    saveas(gcf,['Sticky_WiSE_05m_rho_Mesh_01_a_01.fig'])

    figure
    fig=gcf;
%     fig.Position(3) = 400;
%     fig.Position(4) = 400;
    contourf(Phi,DPhi,phi(:,:,1)',250,'edgecolor','none')
    %plot(Phi,phi(:,:,3)','linewidth',1.5)
    hold on
    %colorbar
    %title('$\bar{\phi}_0$', 'interpreter','latex', 'fontsize', 16)
    %title('$\phi_0$ dependence on $\tilde{\nabla}u$ \& $u$', 'interpreter','latex', 'fontsize', 16)
    ylabel('$\tilde{\nabla}u$', 'interpreter','latex', 'fontsize', 16)
    xlabel('$u$', 'interpreter','latex', 'fontsize', 16)
    cb=colorbar();
    ylabel(cb,'$\bar{\phi}_+$', 'interpreter','latex', 'fontsize', 16)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultTextInterpreter', 'latex')
    ax = gca;
    ax.FontSize = 16;
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'Units','normalized', 'Position',[0.1 0.15 0.675 0.8]);
    set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','on');
    exportgraphics(gcf,['Sticky_WiSE_05m_phic_Mesh_01_a_01.eps'],'Resolution',600)
    exportgraphics(gcf,['Sticky_WiSE_05m_phic_Mesh_01_a_01.jpeg'],'Resolution',600)
    saveas(gcf,['Sticky_WiSE_05m_phic_Mesh_01_a_01.fig'])

    figure 
    fig=gcf;
%     fig.Position(3) = 400;
%     fig.Position(4) = 400;
    contourf(Phi,DPhi,phi(:,:,2)',250,'edgecolor','none')
    %plot(Phi,phi(:,:,3)','linewidth',1.5)
    hold on
    %colorbar
    %title('$\bar{\phi}_0$', 'interpreter','latex', 'fontsize', 16)
    %title('$\phi_0$ dependence on $\tilde{\nabla}u$ \& $u$', 'interpreter','latex', 'fontsize', 16)
    ylabel('$\tilde{\nabla}u$', 'interpreter','latex', 'fontsize', 16)
    xlabel('$u$', 'interpreter','latex', 'fontsize', 16)
    cb=colorbar();
    ylabel(cb,'$\bar{\phi}_-$', 'interpreter','latex', 'fontsize', 16)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultTextInterpreter', 'latex')
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'Units','normalized', 'Position',[0.1 0.15 0.675 0.8]);
    set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','on');
    exportgraphics(gcf,['Sticky_WiSE_05m_phia_Mesh_01_a_01.eps'],'Resolution',600)
    exportgraphics(gcf,['Sticky_WiSE_05m_phia_Mesh_01_a_01.jpeg'],'Resolution',600)
    saveas(gcf,['Sticky_WiSE_05m_phia_Mesh_01_a_01.fig'])
    % saveas(gcf,'Sticky_WiSE_05m_phia_Mesh','jpeg')
    % saveas(gcf,'Sticky_WiSE_05m_phia_Mesh','epsc')
    % saveas(gcf,'Sticky_WiSE_05m_phia_Mesh','fig')

    figure
    fig=gcf;
%     fig.Position(3) = 400;
%     fig.Position(4) = 400;
    contourf(Phi,DPhi,phi(:,:,3)',250,'edgecolor','none')
    %plot(Phi,phi(:,:,3)','linewidth',1.5)
    hold on
    %colorbar
    %title('$\bar{\phi}_0$', 'interpreter','latex', 'fontsize', 16)
    %title('$\phi_0$ dependence on $\tilde{\nabla}u$ \& $u$', 'interpreter','latex', 'fontsize', 16)
    ylabel('$\tilde{\nabla}u$', 'interpreter','latex', 'fontsize', 16)
    xlabel('$u$', 'interpreter','latex', 'fontsize', 16)
    cb=colorbar();
    ylabel(cb,'$\bar{\phi}_0$', 'interpreter','latex', 'fontsize', 16)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultTextInterpreter', 'latex')
    ax = gca;
    ax.FontSize = 16;
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'Units','normalized', 'Position',[0.1 0.15 0.675 0.8]);
    set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','on');
    exportgraphics(gcf,['Sticky_WiSE_05m_phi0_Mesh_01_a_01.eps'],'Resolution',600)
    exportgraphics(gcf,['Sticky_WiSE_05m_phi0_Mesh_01_a_01.jpeg'],'Resolution',600)
    saveas(gcf,['Sticky_WiSE_05m_phi0_Mesh_01_a_01.fig'])
    % saveas(gcf,'Sticky_WiSE_05m_phi0_Mesh','jpeg')
    % saveas(gcf,'Sticky_WiSE_05m_phi0_Mesh','epsc')
    % saveas(gcf,'Sticky_WiSE_05m_phi0_Mesh','fig')

    figure
    fig=gcf;
%     fig.Position(3) = 400;
%     fig.Position(4) = 400;
    contourf(Phi,DPhi,(eta_a)',250,'edgecolor','none')
    %plot(Phi,phi(:,:,3)','linewidth',1.5)
    hold on
    % colorbar
    %title('$\bar{p}_{+0}$', 'interpreter','latex', 'fontsize', 16)
    %title('$\phi_0$ dependence on $\tilde{\nabla}u$ \& $u$', 'interpreter','latex', 'fontsize', 16)
    ylabel('$\tilde{\nabla}u$', 'interpreter','latex', 'fontsize', 16)
    xlabel('$u$', 'interpreter','latex', 'fontsize', 16)
    cb=colorbar();
    ylabel(cb,'$\bar{p}_{+0}$', 'interpreter','latex', 'fontsize', 16)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultTextInterpreter', 'latex')
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'Units','normalized', 'Position',[0.1 0.15 0.675 0.8]);
    set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','on');
    exportgraphics(gcf,['Sticky_WiSE_05m_pp0_Mesh_01_a_01.eps'],'Resolution',600)
    exportgraphics(gcf,['Sticky_WiSE_05m_pp0_Mesh_01_a_01.jpeg'],'Resolution',600)
    saveas(gcf,['Sticky_WiSE_05m_pp0_Mesh_01_a_01.fig'])
    % saveas(gcf,'Sticky_WiSE_05m_pp0_Mesh','jpeg')
    % saveas(gcf,'Sticky_WiSE_05m_pp0_Mesh','epsc')
    % saveas(gcf,'Sticky_WiSE_05m_pp0_Mesh','fig')

    figure
    fig=gcf;
%     fig.Position(3) = 400;
%     fig.Position(4) = 400;
    contourf(Phi,DPhi,(1-eta_a)',250,'edgecolor','none')
    %plot(Phi,phi(:,:,3)','linewidth',1.5)
    hold on
    % colorbar
    %title('$\bar{p}_{+0}$', 'interpreter','latex', 'fontsize', 16)
    %title('$\phi_0$ dependence on $\tilde{\nabla}u$ \& $u$', 'interpreter','latex', 'fontsize', 16)
    ylabel('$\tilde{\nabla}u$', 'interpreter','latex', 'fontsize', 16)
    xlabel('$u$', 'interpreter','latex', 'fontsize', 16)
    cb=colorbar();
    ylabel(cb,'$\bar{p}_{+-}$', 'interpreter','latex', 'fontsize', 16)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultTextInterpreter', 'latex')
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'Units','normalized', 'Position',[0.1 0.15 0.675 0.8]);
    set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','on');
    exportgraphics(gcf,['Sticky_WiSE_05m_ppm_Mesh_01_a_01.eps'],'Resolution',600)
    exportgraphics(gcf,['Sticky_WiSE_05m_ppm_Mesh_01_a_01.jpeg'],'Resolution',600)
    saveas(gcf,['Sticky_WiSE_05m_ppm_Mesh_01_a_01.fig'])
    % saveas(gcf,'Sticky_WiSE_05m_ppm_Mesh','jpeg')
    % saveas(gcf,'Sticky_WiSE_05m_ppm_Mesh','epsc')
    % saveas(gcf,'Sticky_WiSE_05m_ppm_Mesh','fig')

    figure
    fig=gcf;
%     fig.Position(3) = 400;
%     fig.Position(4) = 400;
    contourf(Phi,DPhi,(pa_a)',250,'edgecolor','none')
    %plot(Phi,phi(:,:,3)','linewidth',1.5)
    hold on
    % colorbar
    %title('$\bar{p}_{+0}$', 'interpreter','latex', 'fontsize', 16)
    %title('$\phi_0$ dependence on $\tilde{\nabla}u$ \& $u$', 'interpreter','latex', 'fontsize', 16)
    ylabel('$\tilde{\nabla}u$', 'interpreter','latex', 'fontsize', 16)
    xlabel('$u$', 'interpreter','latex', 'fontsize', 16)
    cb=colorbar();
    ylabel(cb,'$\bar{p}_{-+}$', 'interpreter','latex', 'fontsize', 16)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultTextInterpreter', 'latex')
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'Units','normalized', 'Position',[0.1 0.15 0.675 0.8]);
    set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','on');
    exportgraphics(gcf,['Sticky_WiSE_05m_pmp_Mesh_01_a_01.eps'],'Resolution',600)
    exportgraphics(gcf,['Sticky_WiSE_05m_pmp_Mesh_01_a_01.jpeg'],'Resolution',600)
    saveas(gcf,['Sticky_WiSE_05m_pmp_Mesh_01_a_01.fig'])
    % saveas(gcf,'Sticky_WiSE_05m_pmp_Mesh','jpeg')
    % saveas(gcf,'Sticky_WiSE_05m_pmp_Mesh','epsc')
    % saveas(gcf,'Sticky_WiSE_05m_pmp_Mesh','fig')

    figure
    fig=gcf;
%     fig.Position(3) = 400;
%     fig.Position(4) = 400;
    contourf(Phi,DPhi,(p0_a)',250,'edgecolor','none')
    %plot(Phi,phi(:,:,3)','linewidth',1.5)
    hold on
    % colorbar
    %title('$\bar{p}_{+0}$', 'interpreter','latex', 'fontsize', 16)
    %title('$\phi_0$ dependence on $\tilde{\nabla}u$ \& $u$', 'interpreter','latex', 'fontsize', 16)
    ylabel('$\tilde{\nabla}u$', 'interpreter','latex', 'fontsize', 16)
    xlabel('$u$', 'interpreter','latex', 'fontsize', 16)
    cb=colorbar();
    ylabel(cb,'$\bar{p}_{0+}$', 'interpreter','latex', 'fontsize', 16)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultTextInterpreter', 'latex')
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'Units','normalized', 'Position',[0.1 0.15 0.675 0.8]);
    set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','on');
    exportgraphics(gcf,['Sticky_WiSE_05m_p0p_Mesh_01_a_01.eps'],'Resolution',600)
    exportgraphics(gcf,['Sticky_WiSE_05m_p0p_Mesh_01_a_01.jpeg'],'Resolution',600)
    saveas(gcf,['Sticky_WiSE_05m_p0p_Mesh_01_a_01.fig'])
    % saveas(gcf,'Sticky_WiSE_05m_p0p_Mesh','jpeg')
    % saveas(gcf,'Sticky_WiSE_05m_p0p_Mesh','epsc')
    % saveas(gcf,'Sticky_WiSE_05m_p0p_Mesh','fig')


    figure
    fig=gcf;
%     fig.Position(3) = 400;
%     fig.Position(4) = 400;
    contourf(Phi,DPhi,tau_a',250,'edgecolor','none')
    %plot(Phi,phi(:,:,3)','linewidth',1.5)
    hold on
    %colorbar
    %title('$\tau$', 'interpreter','latex', 'fontsize', 16)
    %title('$\phi_0$ dependence on $\tilde{\nabla}u$ \& $u$', 'interpreter','latex', 'fontsize', 16)
    ylabel('$\tilde{\nabla}u$', 'interpreter','latex', 'fontsize', 16)
    xlabel('$u$', 'interpreter','latex', 'fontsize', 16)
    cb=colorbar();
    ylabel(cb,'$\tau$', 'interpreter','latex', 'fontsize', 16)%,'
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultTextInterpreter', 'latex')
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'Units','normalized', 'Position',[0.1 0.15 0.675 0.8]);
    set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','on');
    exportgraphics(gcf,['Sticky_WiSE_05m_tau_Mesh_01_a_01.eps'],'Resolution',600)
    exportgraphics(gcf,['Sticky_WiSE_05m_tau_Mesh_01_a_01.jpeg'],'Resolution',600)
    saveas(gcf,['Sticky_WiSE_05m_tau_Mesh_01_a_01.fig'])
    % saveas(gcf,'Sticky_WiSE_05m_tau_Mesh','jpeg')
    % saveas(gcf,'Sticky_WiSE_05m_tau_Mesh','epsc')
    % saveas(gcf,'Sticky_WiSE_05m_tau_Mesh','fig')

    figure
    fig=gcf;
    %fig.Position(3) = 400;
    %fig.Position(4) = 400;
    contourf(Phi,DPhi,epsr_a',250,'edgecolor','none')
    %plot(Phi,phi(:,:,3)','linewidth',1.5)
    hold on
    %colorbar
    %title('$\tilde{\epsilon}$', 'interpreter','latex', 'fontsize', 16)
    %title('$\phi_0$ dependence on $\tilde{\nabla}u$ \& $u$', 'interpreter','latex', 'fontsize', 16)
    ylabel('$\tilde{\nabla}u$', 'interpreter','latex', 'fontsize', 16)
    xlabel('$u$', 'interpreter','latex', 'fontsize', 16)
    cb=colorbar();
    ylabel(cb,'$\tilde{\epsilon}$', 'interpreter','latex', 'fontsize', 16)%,'Rotation',270)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultTextInterpreter', 'latex')
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'Units','normalized', 'Position',[0.1 0.15 0.675 0.8]);
    set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','on');
    exportgraphics(gcf,['Sticky_WiSE_05m_Nond_DieCon_Mesh_01_a_01.eps'],'Resolution',600)
    exportgraphics(gcf,['Sticky_WiSE_05m_Nond_DieCon_Mesh_01_a_01.jpeg'],'Resolution',600)
    saveas(gcf,['Sticky_WiSE_05m_Nond_DieCon_Mesh_01_a_01.fig'])
    % saveas(gcf,'Sticky_WiSE_05m_Nond_DieCon_Mesh','jpeg')
    % saveas(gcf,'Sticky_WiSE_05m_Nond_DieCon_Mesh','epsc')
    % saveas(gcf,'Sticky_WiSE_05m_Nond_DieCon_Mesh','fig')

    figure
    fig=gcf;
%     fig.Position(3) = 400;
%     fig.Position(4) = 400;
    contourf(Phi,DPhi,prox_gel',250,'edgecolor','none')
    hold on
    %colorbar
    %title('$\tilde{\epsilon}$', 'interpreter','latex', 'fontsize', 16)
    %title('$\phi_0$ dependence on $\tilde{\nabla}u$ \& $u$', 'interpreter','latex', 'fontsize', 16)
    ylabel('$\tilde{\nabla}u$', 'interpreter','latex', 'fontsize', 16)
    xlabel('$u$', 'interpreter','latex', 'fontsize', 16)
    cb=colorbar();
    ylabel(cb,'Proximity to Gelation', 'interpreter','latex', 'fontsize', 16)%,'Rotation',270)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultTextInterpreter', 'latex')
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'Units','normalized', 'Position',[0.1 0.15 0.675 0.8]);
    set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','on');
    exportgraphics(gcf,['Sticky_WiSE_05m_Proximity_to_Gelation_Mesh_01_a_01.eps'],'Resolution',600)
    exportgraphics(gcf,['Sticky_WiSE_05m_Proximity_to_Gelation_Mesh_01_a_01.jpeg'],'Resolution',600)
    saveas(gcf,['Sticky_WiSE_05m_Proximity_to_Gelation_Mesh_01_a_01.fig'])
    % saveas(gcf,'Sticky_WiSE_05m_Proximity_to_Gelation_Mesh','jpeg')
    % saveas(gcf,'Sticky_WiSE_05m_Proximity_to_Gelation_Mesh','epsc')
    % saveas(gcf,'Sticky_WiSE_05m_Proximity_to_Gelation_Mesh','fig')

end

function [rho,ip_vf,agg_vf,sum_cl,sum_cm,sum_cs] = WiSE_sti_theory_func(fa,fc,Lbar,Xp,Yp,zp,EP)

cl = zeros([EP EP]);
cm = zeros([EP EP]);
cs = zeros([EP EP]);

for J = 0:EP
    for K = 0:EP
        l = J+K;
        %rho = J-K;
        %q = J-K;
        if l > 0 && fc * J - J - K + 1 >= 0 && fa * K - K - J + 1 >=0
            Temp = zp*stockmeyer_prob_4sum(J,K,Xp,Yp,fc,fa);

            cl(J+1,K+1) = J*Temp;%/2;
            cm(J+1,K+1) = K*Temp;%/2;
            cs(J+1,K+1) = (fc*J-J-K+1)*Temp;

        elseif J == 0 && K == 0
            cs(J+1,K+1) = zp*Lbar;
        end
    end
    
end
ip_vf = cl(2,2)+cm(2,2);
sum_cl = sum(sum(cl(:,:),1),2);
sum_cm = sum(sum(cm(:,:),1),2);
sum_cs = sum(sum(cs(:,:),1),2);
agg_vf = sum_cl + sum_cm - ip_vf - cl(2,1) - cm(1,2);
rho = sum_cl- sum_cm;
end

function SN = stockmeyer_prob_4sum(l,m,XX,YY,fc,fa)

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

function B = sinch(x)
    if abs(x) > 1e-3
        B = sinh(x)/x;
    else
        B = 1 + x^2 / 6 + x^4 / 120 + x^6 / 5040;
    end
end

function pc = f_pc(phi_0,psi_p,psi_m,L)
    pc = (phi_0 - psi_p + L*(psi_p + psi_m) - sqrt(4*phi_0*psi_p*(L-1) + (L*(psi_m - psi_p) + psi_p + phi_0)^2))/(2*(L-1)*psi_p);
end

function pa = f_pa(phi_0,psi_p,psi_m,L)
    pa = (phi_0 - psi_p + L*(psi_p + psi_m) - sqrt(4*phi_0*psi_p*(L-1) + (L*(psi_m - psi_p) + psi_p + phi_0)^2))/(2*(L-1)*psi_m);
end

function pc0 = f_pc0(phi_0,psi_p,psi_m,L)
    pc0 = (phi_0 + psi_p + L*(-psi_p + psi_m) - sqrt(4*phi_0*psi_p*(L-1) + (L*(psi_m - psi_p) + psi_p + phi_0)^2))/(2*(1-L)*psi_p);
end

function p0 = f_p0(phi_0,psi_p,psi_m,L)
    p0 = (phi_0 + psi_p + L*(-psi_p + psi_m) - sqrt(4*phi_0*psi_p*(L-1) + (L*(psi_m - psi_p) + psi_p + phi_0)^2))/(2*(1-L)*phi_0);
end

function L_f = f_L_f(p,up)
    %Make loop for vectorized up
    L_f=up;
    for i=1:length(up)
        if abs(up(i)) > 1e-3
            L_f(i) = (coth(p*up(i)) - 1/( p*up(i) ))/(up(i));
        else
            k = p * up(i);
            L_f(i) = p * (1/3 - k^2/45 + 2*k^4/945 - k^6/4725); 
        end
    end
end

function t_epsr = f_t_epsr(cpb,cmb,u,du,p,phi_0,p0p,select)
    t_epsr = cpb + cmb + phi_0.*(1-p0p).* p * f_L_f(p,du);
end
