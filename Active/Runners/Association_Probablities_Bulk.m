%Probability Solver
close all
clear all

%Physical Constants
T = 300; 
Na = 6.022e23; 
kB = 1.38e-23;  
e = 1.6e-19;
beta = 1/(kB*T);

%General Systems Properties and Bulk Quantities
fa = 3;
fc = 4;
xic = 0.4;
xia = 10.8;
MMw=18.02/1000; %kg/mol molar mass of water
L = 0.2527; %Lpm/Lpo; %Ion-Association Constant  

mol_f_ref = 1/((1/MMw/0.5)+1);
mol_f = [mol_f_ref:0.001:0.02-0.001,0.02:0.0025:.35];
molality = 1/MMw.*(1./(1./mol_f-1));
for i = 1:length(mol_f)
    phi_p_b(i) = xic/(1/mol_f(i)-1+xic+xia);
    phi_m_b(i) = xia/xic*phi_p_b(i);
    phi_0_b(i) = 1 - phi_p_b(i) - phi_m_b(i);
        
    pc(i) = f_pc(phi_0_b(i),fc*phi_p_b(i)/xic,fa*phi_m_b(i)/xia,L);
    pa(i) = f_pa(phi_0_b(i),fc*phi_p_b(i)/xic,fa*phi_m_b(i)/xia,L);
    pc0(i) = f_pc0(phi_0_b(i),fc*phi_p_b(i)/xic,fa*phi_m_b(i)/xia,L);
    p0(i) = f_p0(phi_0_b(i),fc*phi_p_b(i)/xic,fa*phi_m_b(i)/xia,L);

end

%Bulk Association Probabilities:
figure()
plot(molality,pc,'r-','linewidth',1.5,'DisplayName','$p_{+-}$');
hold on
plot(molality,pc0,'r--','linewidth',1.5,'DisplayName','$p_{+0}$');
plot(molality,pa,'b:','linewidth',1.5,'DisplayName','$p_{-+}$');
plot(molality,p0,'k-.','linewidth',1.5,'DisplayName','$p_{0+}$');
xlim([0,21])
ylim([-0.025,1.175])
box on
xlabel(['LiTFSI molality, m'],'interpreter','latex', 'fontsize', 16)
ylabel('Association probability, $p_{ij}$','interpreter','latex', 'fontsize', 16)
pbaspect([(1+sqrt(5))/2,1,1])
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'DefaultTextInterpreter', 'latex');
leg=legend('location','north','Orientation','horizontal','Interpreter','latex','Color','none','Box','off');
leg.ItemTokenSize = [18,18];
yticks(0:.2:1);
g=gcf;
g.Renderer='painters';
set(gca, 'Xcolor', 'k');
set(gca, 'Ycolor', 'k');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1.5);
set(gca, 'Layer', 'Top');
set(gca, 'color', 'none');
set(g,'InvertHardcopy','on');
exportgraphics(gcf,['Sticky_WiSE_Bulk_21m_L_assoc_v1.eps'],'Resolution',600)
exportgraphics(gcf,['Sticky_WiSE_Bulk_21m_L_assoc_v1.jpeg'],'Resolution',600)
saveas(gcf,['Sticky_WiSE_Bulk_21m_L_assoc_v1.fig'])
    
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