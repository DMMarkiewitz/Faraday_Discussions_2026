function Implicit_Poisson_Solver
    clear all
    close all
    
    %Key File Inputs
    mm = 21; %Molality
    fig_save_name = 'WiSE_21m_fp4_sticky_cl_dist_bulk';%Figure Save name

    %Physical constants
    T = 300; 
    Na = 6.022e23; 
    kB = 1.38e-23;  
    e = 1.6e-19;
    beta = 1/(kB*T);
    EP = 50;
    
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

    for i_c = 1
        pc = f_pc(phi_0_b,psi_p_b,psi_m_b,L);
        pa = f_pa(phi_0_b,psi_p_b,psi_m_b,L);
        pc0 = f_pc0(phi_0_b,psi_p_b,psi_m_b,L);
        p0 = f_p0(phi_0_b,psi_p_b,psi_m_b,L);
        Lbar = L;
        zp = phi_0_b * (1-p0) / L;
        Xp = fc/xic * phi_p_b * (1-pc)^fc / zp;
        Yp = fa/xia * phi_m_b * (1-pa)^fa / zp;
        [clms,cubeAggLsdv0] = WiSE_sti_theory_func_count(fa,fc,xic,xia,L,Xp,Yp,zp,EP);
        AggLs(i_c) = (cubeAggLsdv0*v0)^(1/3);
        clms_a{i_c} = clms;
    end

    figure()
    temp=log10(clms_a{1});
    temp(temp<-8)=nan;
    fig=pcolor(-0.5:size(temp,1)-1.5,-0.5:size(temp,2)-1.5,temp');
    set(fig, 'EdgeColor', 'none')
    clim([-8,0])
    xlim([0,20])
    ylim([0,20])
    cb = colorbar("LineWidth",1.5);
    ylabel(cb,'Log$_{10}$($c_{lms}$)', 'interpreter','latex', 'fontsize', 16)
    hold on
    plot([0,50],[0,50],'--k','LineWidth',2)
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

    exportgraphics(g,[fig_save_name,'.eps'],'Resolution',600)
    exportgraphics(g,[fig_save_name,'.jpeg'],'Resolution',600)
    saveas(g,[fig_save_name,'.fig'])

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