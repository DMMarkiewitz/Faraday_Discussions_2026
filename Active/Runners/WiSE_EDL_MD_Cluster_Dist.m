function WiSE_EDL_esorption
%One needs to manuelly change the plot's save name when analysing
%individual surface charge densities and molalities
clear all
close all
salt = 'LiTFSI';
qs = 0.2;%[0,0.01,0.05,0.1,0.15,0.2]; %[0,0.05,0.1,0.15,0.2];
mm = 21; %molality %10;
Lx = 32.624;
Ly = 32.624;
Lz = 265.248;
L = [Lx,Ly,Lz];
si = sqrt(10);
gw = 1;%if gw = 1 we are not smoothing the data
gw_p = 1;%if gw_p = 1 we are not smoothing the data
nbin = 650; %Number of bins

%Distance for spatial averaging
sp_A = 5; %In angstroms

if isequal(salt,'LiTFSI')
    nc = 1;
    na = 15;
    nao = 4;
elseif isequal(salt,'LiOTF')
    nc = 1;
    na = 8;
    nao = 3;
elseif isequal(salt,'LiFSI')
    nc = 1;
    na = 9;
    nao = 4;
end

nwo = 1;

sep = 2.7; % separatrix for participating/spectating ions
cut_off = 2.7; % cut off for coordination

%association numbers
fp = 4;
fm = 3;
fs = 1;

%Ion association constants:
Lpm_p = 0;
Lpm_m = 0;
Lp0_p = 0;
Lp0_0 = 0;
for i = 1 : length(qs)
    filePath = ['EDL_',num2str(qs(i)),'.xyz'];
    % [csalt,data]=readDumpFile(filePath,1,salt,mm,sep,cut_off,L);
    savename = ['S_EDL_paper26_','sep_',num2str(sep),'_',num2str(qs(i)),'.mat'];
    %save(savename,'data','-v7.3')
    load(savename,'data')
    csalt = getconc(filePath,salt,mm,L);

    csolvent = 55.556/mm*csalt;
    Label = [data.label];
    Labels = [data.label_s];
    Label_anion = [data.label_a];
    
    %Reduces the dimensionality of the problem to solely z axis works well as
    %our system is nice i.e. the x and y boundaries do not vary with zcut_off
    Cat_zpos = [data.com_c];
    An_zpos = [data.com_a];%center of mass of whole anion molecule
    Sol_zpos = [data.com_s];
    
    %Defining bounds
    Lzl = min([min(min(Cat_zpos)),min(min(An_zpos)),min(min(Sol_zpos))])*0.9999;
    Lzu = max([max(max(Cat_zpos)),max(max(An_zpos)),max(max(Sol_zpos))])*(1+1e-4);
    ed = unique([flip(Lzl:-(Lzu-Lzl)/(nbin-1):0) linspace(Lzl,Lzu,nbin) (Lzu:(Lzu-Lzl)/(nbin-1):Lz)]);
    x = ed(2:end) -(ed(2)-ed(1))/2; %Finds the center nodes
    %Defining the active indexes for the interior data
    active_index = length(Lzl:-(Lzu-Lzl)/(nbin-1):0):length(unique([flip(Lzl:-(Lzu-Lzl)/(nbin-1):0) linspace(Lzl,Lzu,nbin)]));
    
    %Parameters for conversion to volume fractions and getting plots with
    %distance in terms of kappa x
    if isequal(salt,'LiTFSI')
        xip = 0.4;
        xim = 10.8;
        omeg = (xip+xim)*length(An_zpos)+length(Sol_zpos);
        Vbox = Lx*Ly*(Lzu-Lzl)*1e-30;
        v0 = Vbox/omeg; %m^3 volume backcalculated from vbox and Omega,
        T = 300; 
        NA = 6.022e23; 
        kB = 1.38e-23;  
        e = 1.6e-19;
        beta = 1/(kB*T);
        eps_r = 10.1; %So our debye length used to non-dimensionalize is same as our theory
        eps=eps_r*8.85e-12;
        phi_0_b = length(Sol_zpos)/omeg;
        phi_p_b = (1-phi_0_b)/(1+xim/xip);
        phi_m_b = 1-phi_p_b-phi_0_b;
        AA2kappa = 1/(sqrt(v0*eps/(beta*e^2*(phi_p_b/xip + phi_m_b/xim)))*1e10);%Angstroms to kappa
    elseif isequal(salt,'LiOTF')
        error('Not yet implemented')
    elseif isequal(salt,'LiFSI')
        error('Not yet implemented')
    end

    %Plotting the coordination of Li with z location
    %As there is values associated with the positions we need to first
    %parse our data and do our summing ourself
    Li_coord_an = zeros(length(x),1); %Preinitialize Li-array, coords are the x bins coords
    Li_coord_w = zeros(length(x),1); %Preinitialize Li-array, coords are the x bins coords
    n_length = length(data); %How many time steps we are averaging our data over
    mean_Li_c_an = 0;
    mean_Li_c_w = 0;
    Li_c_an_s = @(x,sigma) 0;
    Li_c_w_s = @(x,sigma) 0;
    for n = 1:n_length
        Li_c_an_s = @(x,sigma) Li_c_an_s(x,sigma) + sum(data(n).coord_li_an(:,2) .* exp(-((x-data(n).coord_li_an(:,1))/sigma).^2 / 2) / (sigma*sqrt(2*pi)))/n_length;
        Li_c_w_s = @(x,sigma) Li_c_w_s(x,sigma) + sum(data(n).coord_li_w(:,2) .* exp(-((x-data(n).coord_li_w(:,1))/sigma).^2 / 2) / (sigma*sqrt(2*pi)))/n_length;
        for ii = 1:length(x) % for loop running over the bins
            inx1 = (ed(ii) <= data(n).coord_li_an(:,1)) & (data(n).coord_li_an(:,1) <= ed(ii+1));
            inx2 = (ed(ii) <= data(n).coord_li_w(:,1)) & (data(n).coord_li_w(:,1) <= ed(ii+1));
            Li_coord_an(ii) = Li_coord_an(ii) + sum(data(n).coord_li_an(inx1,2))/(n_length*max(1,sum(inx1))); %Summing coordination number of Li in this bin and weigthing it unifromly amond the timesteps and how many species are present
            Li_coord_w(ii) = Li_coord_w(ii) + sum(data(n).coord_li_w(inx2,2))/(n_length*max(1,sum(inx2))); %Summing coordination number of Li in this bin and weigthing it unifromly amond the timesteps steps and how many species are present
        end
        mean_Li_c_an = mean_Li_c_an + mean(data(n).coord_li_an(:,2))/n_length;
        mean_Li_c_w = mean_Li_c_w + mean(data(n).coord_li_w(:,2))/n_length;
    end
    %We want to normalize it so the value within the electrode is accurate
    %for the global average
    Li_c_an_s = @(x,sigma) mean_Li_c_an*Li_c_an_s(x,sigma)/(integral(@(x) Li_c_an_s(x,sigma),Lzl,Lzu)/(Lzu-Lzl));
    Li_c_w_s = @(x,sigma) mean_Li_c_w*Li_c_w_s(x,sigma)/(integral(@(x) Li_c_w_s(x,sigma),Lzl,Lzu)/(Lzu-Lzl));

    %Calculations for association probabilities  
    %As there is values associated with the positions we need to first
    %parse our data and do our summing ourself
    Op = zeros(length(x),1); %Preinitialize Op-array, coords are the x bins coords
    nc = Op; %Preinitialize nc-array, coords are the x bins coords
    na = Op; %Preinitialize na-array, coords are the x bins coords
    ns = Op; %Preinitialize ns-array, coords are the x bins coords
    ns_5 = zeros(2,1);
    nsb_5 = zeros(2,1);
    nsb = Op; %Preinitialize nsb-array, coords are the x bins coords
    n_length = length(data); %How many time steps we are averaging our data over
    n_sol = Op; %Tracks how many timesteps the box is truely active for i.e. have one particle of their class present in
    for n = 1:n_length
        Ncat = histcounts(data(n).com_c,ed);
        Nan = histcounts(data(n).com_a,ed);
        Nsol = histcounts(data(n).xyz_ow(:,3),ed);
        nc = nc + Ncat'/n_length;
        na = na + Nan'/n_length;
        ns = ns + Nsol'/n_length;
        for ii = 1:length(x) % for loop running over the bins
            inx2 = (ed(ii) <= data(n).Opassoc(:,1)) & (data(n).Opassoc(:,1) <= ed(ii+1));   
            nsb(ii) = nsb(ii) + sum(data(n).Opassoc(inx2,2)>0)/n_length; %Average number of bound water molecules, >0 converts all bound species into weight 1 even if they are twice bound as that is current methodology        
        end
        lzl_a = (Lzl <= data(n).Opassoc(:,1)) & (data(n).Opassoc(:,1) <= Lzl+sp_A);
        lzu_a = (Lzu-sp_A <= data(n).Opassoc(:,1)) & (data(n).Opassoc(:,1) <= Lzu);
        nsb_5(1) = nsb_5(1) + sum(data(n).Opassoc(lzl_a,2)>0)/n_length;
        nsb_5(2) = nsb_5(2) + sum(data(n).Opassoc(lzu_a,2)>0)/n_length;
        lzl_a = (Lzl <= data(n).xyz_ow(:,3)) & (data(n).xyz_ow(:,3) <= Lzl+sp_A);
        lzu_a = (Lzu-sp_A <= data(n).xyz_ow(:,3)) & (data(n).xyz_ow(:,3) <= Lzu);
        ns_5(1) = ns_5(1) + sum(data(n).xyz_ow(lzl_a,3)>0)/n_length;
        ns_5(2) = ns_5(2) + sum(data(n).xyz_ow(lzu_a,3)>0)/n_length;
    end
    nsf = ns - nsb; %Extracting the free molecules
    nsf_5 = ns_5 - nsb_5; %Extracting the free molecules in a 5 A box from surface

    %Calculating the pij(x) temporally averaged
    Op = Op./(fs*n_sol);

    %Determining cutting points
    [xmin,ix] = min([min(x(nc~=0)),min(x(na~=0)),min(x(ns~=0))]);
    [xmax,ixma] = max([max(x(nc~=0)),max(x(na~=0)),max(x(ns~=0))]);

    if ix==1
        act_mi = nc~=0;
    elseif ix==2
        act_mi = na~=0;
    else
        act_mi = ns~=0;
    end

    if ixma==1
        act_ma = nc~=0;
    elseif ixma==2
        act_ma = na~=0;
    else
        act_ma = ns~=0;
    end
    
    act_i = or(act_mi,act_ma);

    xc=x(act_i);

    %Finding the region of expected model breakdown
    [gl,~] = max([min(x(nc~=0)),min(x(na~=0)),min(x(ns~=0))]);
    [gu,~] = min([max(x(nc~=0)),max(x(na~=0)),max(x(ns~=0))]);
    xminvf = xmin;
    xmaxvf = xmax;

    %Points of interest
    [~,xi10p]= min(abs((x-xmin)*AA2kappa - 10));
    [~,xi5p]= min(abs((x-xmin)*AA2kappa - 5));
    [~,xi10m]= min(abs((x(act_i)-xmin)*AA2kappa - 10));
    [~,xi5m]= min(abs((x(act_i)-xmin)*AA2kappa - 5));
    
    %Making a list of all the types of clusters
    clusterlist = string(data(1).typesofcluster); 
    for ia = 2:n
       clusterlist = cat(1,clusterlist,string(data(ia).typesofcluster));
    end
    
    %Counting the occurances of the clusters
    Count = zeros(length(Cat_zpos)+1,length(Cat_zpos)+1,length(x));
    Count_5 = zeros(length(Cat_zpos)+1,length(Cat_zpos)+1,2);
    max_l = 0;
    max_m = 0;
    for inx = 1:n_length
        if max(data(inx).cluster(:,3)) > max_l
            max_l = max(data(inx).cluster(:,3));
        end
        if max(data(inx).cluster(:,4)) > max_m
            max_m = max(data(inx).cluster(:,4));
        end
        for ii = 1:length(x) % for loop running over the bins
            inx1 = (ed(ii) <= data(inx).cluster(:,2)) & (data(inx).cluster(:,2) <= ed(ii+1));
            Count(:,:,ii) = Count(:,:,ii) + sparse(data(inx).cluster(inx1,3)+1,data(inx).cluster(inx1,4)+1,1./(data(inx).cluster(inx1,3)+data(inx).cluster(inx1,4)),length(Cat_zpos)+1,length(Cat_zpos)+1)/length(data);
        end
        lzl_a = (Lzl <= data(inx).cluster(:,2)) & (data(inx).cluster(:,2) <= Lzl+sp_A);
        lzu_a = (Lzu-sp_A <= data(inx).cluster(:,2)) & (data(inx).cluster(:,2) <= Lzu);
        Count_5(:,:,1) = Count_5(:,:,1) + sparse(data(inx).cluster(lzl_a,3)+1,data(inx).cluster(lzl_a,4)+1,1./(data(inx).cluster(lzl_a,3)+data(inx).cluster(lzl_a,4)),length(Cat_zpos)+1,length(Cat_zpos)+1)/length(data);
        Count_5(:,:,2) = Count_5(:,:,2) + sparse(data(inx).cluster(lzu_a,3)+1,data(inx).cluster(lzu_a,4)+1,1./(data(inx).cluster(lzu_a,3)+data(inx).cluster(lzu_a,4)),length(Cat_zpos)+1,length(Cat_zpos)+1)/length(data);
    end
    Count=Count*v0/((ed(2)-ed(1))*Lx*Ly*10^-30);
    Count_5=Count_5*v0/(sp_A*Lx*Ly*10^-30);

    %Adding in Count for free solvent as nsf is the average in a bin
    Count(1,1,:)=nsf*v0/((ed(2)-ed(1))*Lx*Ly*10^-30);
    Count_5(1,1,:)=nsf_5*v0/(sp_A*Lx*Ly*10^-30);

    %For the negative side
    Count_m = flip(Count(:,:,act_i),3);

    %Plotting Averaged Bins
    %Positive
    A=log10(Count_5(:,:,1));
    A(A<-8)=nan;
    figure()
    fig=pcolor(-0.5:size(Count_5,1)-1.5,-0.5:size(Count_5,1)-1.5,A')
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
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_p0_2_5A_01.eps','Resolution',600)
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_p0_2_5A_01.jpeg','Resolution',600)
    saveas(g,'WiSE_21m_fp4_sticky_cl_dist_p0_2_5A_01.fig')

    %Negative
    A=log10(Count_5(:,:,2));
    A(A<-8)=nan;
    figure()
    fig=pcolor(-0.5:size(Count_5,1)-1.5,-0.5:size(Count_5,1)-1.5,A')
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
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_m0_2_5A_01.eps','Resolution',600)
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_m0_2_5A_01.jpeg','Resolution',600)
    saveas(g,'WiSE_21m_fp4_sticky_cl_dist_m0_2_5A_01.fig')

    %Plotting Positional Cluster Distributions
    %Postive charge at 5 debye lengths from surface
    A=log10(Count(:,:,xi5p));
    A(A<-8)=nan;
    figure()
    fig=pcolor(-0.5:size(Count,1)-1.5,-0.5:size(Count,1)-1.5,A')
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
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_p0_2_pl_5_01.eps','Resolution',600)
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_p0_2_pl_5_01.jpeg','Resolution',600)
    saveas(g,'WiSE_21m_fp4_sticky_cl_dist_p0_2_pl_5_01.fig')

    %Negative charge at 5 debye lengths from surface
    A=log10(Count_m(:,:,xi5m));
    A(A<-8)=nan;
    figure()
    fig=pcolor(-0.5:size(Count,1)-1.5,-0.5:size(Count,1)-1.5,A')
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
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_m0_2_pl_5_01.eps','Resolution',600)
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_m0_2_pl_5_01.jpeg','Resolution',600)
    saveas(g,'WiSE_21m_fp4_sticky_cl_dist_m0_2_pl_5_01.fig')

    %Positive charge at 10 debye lengths from surface
    A=log10(Count(:,:,xi10p));
    A(A<-8)=nan;
    figure()
    fig=pcolor(-0.5:size(Count,1)-1.5,-0.5:size(Count,1)-1.5,A')
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
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_p0_2_pl_10_01.eps','Resolution',600)
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_p0_2_pl_10_01.jpeg','Resolution',600)
    saveas(g,'WiSE_21m_fp4_sticky_cl_dist_p0_2_pl_10_01.fig')

    %Negative charge at 10 debye lengths from surface
    A=log10(Count_m(:,:,xi10m));
    A(A<-8)=nan;
    figure()
    fig=pcolor(-0.5:size(Count,1)-1.5,-0.5:size(Count,1)-1.5,A')
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
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_m0_2_pl_10_01.eps','Resolution',600)
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_m0_2_pl_10_01.jpeg','Resolution',600)
    saveas(g,'WiSE_21m_fp4_sticky_cl_dist_m0_2_pl_10_01.fig')

    %Plotting Bulk Cluster Distribution
    Count_b = Count(:,:,length(Lzl:-(Lzu-Lzl)/(nbin-1):0):length(unique([flip(Lzl:-(Lzu-Lzl)/(nbin-1):0) linspace(Lzl,Lzu,nbin)])));
    Count_b = mean(Count_b(:,:,ceil(nbin/3):floor(2*nbin/3)),3);

    A=log10(Count_b);
    A(A<-8)=nan;
    figure()
    fig=pcolor(-0.5:size(Count,1)-1.5,-0.5:size(Count,1)-1.5,A')
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
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_b0_2_place_bulk_01.eps','Resolution',600)
    exportgraphics(g,'WiSE_21m_fp4_sticky_cl_dist_b0_2_place_bulk_01.jpeg','Resolution',600)
    saveas(g,'WiSE_21m_fp4_sticky_cl_dist_b0_2_place_bulk_01.fig')

end

end

function csalt=getconc(filePath,salt,mm,L)

fid=fopen(filePath);

N=textscan(fid,'%f',1);
filetext = fileread(filePath); %To accurately detect values will use uniqueness of Li and atoms headers to determine n_c total
n_Li_total = length(regexp(filetext,salt(1:2),'match'));
n_timestep = length(regexp(filetext,'Atoms','match'));

if isequal(salt,'LiTFSI')
    nc = 1;
    na = 15;
    nao = 4;
elseif isequal(salt,'LiOTF')
    nc = 1;
    na = 8;
    nao = 3;
elseif isequal(salt,'LiFSI')
    nc = 1;
    na = 9;
    nao = 4;
end
ns = 3;

n_gold = 4096;
n_tot = N{1}-n_gold;
m = 55.556/mm;
n_ion = n_Li_total/n_timestep; %Number of salt i.e. ionic species is equal to total Li in each simulation
n_solvent = (n_tot - n_ion*(nc+na))/ns;
csalt = n_ion/L(1)/L(2)/L(3)*1e30./6.022e23/1000;

fclose all;
end

function [csalt,data]=readDumpFile(filePath,stepSize,salt,mm,sep,cut_off,L)

fid=fopen(filePath);

N=textscan(fid,'%f',1);
textscan(fid,'%s %s %f',1);
k=1;
filetext = fileread(filePath); %To accurately detect values will use uniqueness of Li and atoms headers to determine n_c total
n_Li_total = length(regexp(filetext,salt(1:2),'match'));
n_timestep = length(regexp(filetext,'Atoms','match'));

if isequal(salt,'LiTFSI')
    nc = 1;
    na = 15;
    nao = 4;
elseif isequal(salt,'LiOTF')
    nc = 1;
    na = 8;
    nao = 3;
elseif isequal(salt,'LiFSI')
    nc = 1;
    na = 9;
    nao = 4;
end
ns = 3; % atoms in solvent molecules, water H2O

n_gold = 4096;
n_tot = N{1}-n_gold;
n_ion = n_Li_total/n_timestep; %Number of salt i.e. ionic species is equal to total Li in each simulation
n_solvent = (n_tot - n_ion*(nc+na))/ns;
csalt = n_ion/L(1)/L(2)/L(3)*1e30./6.022e23/1000;

%load mass from external files, S Bi
mass_all=load(strcat(salt,'.mass'));
mass_cat=sum(mass_all(1:nc));
mass_anion=sum(mass_all(nc+1:end));
cat_massList=mass_all(1:nc);
anion_massList=mass_all(nc+1:end);
tim = {0,0,0};
while tim{3}*10^-6 < 0  %Allowing one to exclude equilibrium period
    textscan(fid,'%s %s %d %f %f %f',n_gold); % grabs all gold atoms
    tmp=textscan(fid,'%s %f %f %f',n_ion*nc); %grabs all cations atoms
    tmp2=textscan(fid,'%s %f %f %f',n_ion*na); %grabs all anion atoms
    tmp3=textscan(fid,'%s %f %f %f',n_solvent*ns); %grabs all water mol
    N=textscan(fid,'%f',1);%Grabs number of atoms in next simulation
    tim = textscan(fid,'%s %s %f',1);
end

while ~isempty(N{1})
    
    for step=1:stepSize
        textscan(fid,'%s %s %d %f %f %f',n_gold); % grabs all gold atoms
        tmp=textscan(fid,'%s %f %f %f',n_ion*nc); %grabs all cations atoms
        tmp2=textscan(fid,'%s %f %f %f',n_ion*na); %grabs all anion atoms
        tmp3=textscan(fid,'%s %f %f %f',n_solvent*ns); %grabs all water mol
        N=textscan(fid,'%f',1);%Grabs number of atoms in next simulation
        textscan(fid,'%s %s %f',1); %Skips over timestep information
        
    end
    % Break statement in case our stepSize is not a mod of the number
    % of timesteps we took thus avoiding errors
    if isempty(N{1}) && length(tmp{1}) == 0
       break 
    end
    
    %Locating center of masses and key motifs
    xyz_cat=horzcat(tmp{2}, tmp{3}, tmp{4});
    xyz_an = horzcat(tmp2{2},tmp2{3},tmp2{4});
    xyz_w = horzcat(tmp3{2},tmp3{3},tmp3{4});
    xyz_oa = xyz_an(strcmp(tmp2{1},'O'),:); %Detects all the association O atoms in the anions
    xyz_ow = xyz_w(strcmp(tmp3{1},'O'),:);  %Detects all the O atoms in the water molecules
    
    %Saving some of the xyz locations to data for potential future usage
    data(k).xyz_cat = xyz_cat;
    data(k).xyz_oa = xyz_oa;
    data(k).xyz_ow = xyz_ow;
    
    % calculate the center-of-mass position of anion, S Bi
    for i=1:n_ion
        indexStart = (i-1)*na+1;
        indexEnd = i*na;
        com_a_x(i) = dot(tmp2{2}(indexStart:indexEnd),anion_massList)/mass_anion;
        com_a_y(i) = dot(tmp2{3}(indexStart:indexEnd),anion_massList)/mass_anion;
        com_a_z(i) = dot(tmp2{4}(indexStart:indexEnd),anion_massList)/mass_anion;
    end
    data(k).com_a = com_a_z';
    data(k).com_c = xyz_cat(:,3);
    data(k).com_s = xyz_ow(:,3);
    
    periodicDist2 = @(x,y)periodicDist(x, y, L, 0);%Z flag is 0 as we only have periodicidity in x and y boundaries of the box

    [~, d] = knnsearch(xyz_oa,xyz_cat,'k',1,'distance',periodicDist2);
    data(k).D = d; %distance between lithium and closest anionic oxygen
    data(k).label = data(k).D < sep;
    
    % Obtain participating/spectating anions
    [~, da] = knnsearch(xyz_cat,xyz_oa,'k',1,'distance',periodicDist2);
    dap = min(reshape(da,[nao,n_ion]),[],1);
    data(k).Da = dap';
    data(k).label_a = data(k).Da < sep;
    
    % Coordination numbers of Spectating Lithium by spectating anion
    [~, dda]=rangesearch(xyz_oa,xyz_cat,sep,'distance',periodicDist2);
    data(k).coorda = cellfun('size',dda,2);
    data(k).type=tmp{1};
    
    % Coordination numbers of Bound and Participating Ions
    [~, dd]=rangesearch(xyz_ow,xyz_cat,cut_off,'distance',periodicDist2);
    data(k).coord = cellfun('size',dd,2);
    data(k).coord_li_w = [xyz_cat(:,3),cellfun('size',dd,2)]; %vertical concatination to make plotting and future work simpler in the future

    % Coordination numbers of Li with anions
    [id, dd]=rangesearch(xyz_oa,xyz_cat,sep,'distance',periodicDist2);
    data(k).coord_li_an = [xyz_cat(:,3),cellfun(@(x) length(unique(ceil(x/4))),id)]; %vertical concatination to make plotting and future work simpler in the future, assuming 1 anion cannot form 2 associations with the cation
    
    %+0 associations
    [~, dd]=rangesearch(xyz_ow,xyz_cat,cut_off,'distance',periodicDist2);
    data(k).pOassoc = [xyz_cat(:,3),cellfun('size',dd,2)]; %vertical concatination to make plotting and future work simpler in the future

    %+- associations
    [id, dd]=rangesearch(xyz_oa,xyz_cat,sep,'distance',periodicDist2);
    data(k).pmassoc = [xyz_cat(:,3),cellfun(@(x) length(unique(ceil(x/4))),id)]; %vertical concatination to make plotting and future work simpler in the future, cannot use simple count as it can bind to multiple O's on one anion and that only counts once

    %0+ associations
    %Done in the simplest manner placing association at oxygen location
    [~, dd]=rangesearch(xyz_cat,xyz_ow,cut_off,'distance',periodicDist2);
    data(k).Opassoc = [xyz_ow(:,3),cellfun('size',dd,2)]; %vertical concatination to make plotting and future work simpler in the future
    
    %-+ associations
    %Done in the simplest manner placing association at com of the anion
    [id, dd]=rangesearch(xyz_cat,xyz_oa,sep,'distance',periodicDist2);
    data(k).mpassoc = zeros(length(data(k).pmassoc),2);
    for ik = 1:length(xyz_cat(:,3))
        data(k).mpassoc(ik,:) = [data(k).com_a(ik),length(unique(cat(2,id{ik*4-3:ik*4})))]; 
    end

    [~, ds]=knnsearch(xyz_cat,xyz_ow,'k',1,'distance',periodicDist2);
    data(k).Ds = ds; % distance between water and closest lithium
    data(k).label_s = data(k).Ds < cut_off; 

    %Currently allowing self-loops
    [idxsa, ~]=rangesearch(xyz_oa,xyz_cat,sep,'distance',periodicDist2);
    p_idxsa = cellfun(@(x) unique(ceil(x/4)),idxsa,'UniformOutput',false); %Removing duplicates
    [idxsw, ~]=rangesearch(xyz_ow,xyz_cat,cut_off,'distance',periodicDist2);
    A = sparse(zeros(2*length(idxsa),2*length(idxsa))); % adjacency matrix, as we have equal amounts of anions and cations
    for q = 1 : length(p_idxsa) %cycle through cations
        for w = 1 : length(p_idxsa{q}) %cycle through paired anions
            ia = p_idxsa{q}(w)+length(idxsa); %This is to say bonding to one oxygen atom in the anion counts as bonding to the anion hence the negative in our mental image is the whole anion not just the oxygen species 
            A(q,ia) = 1;
            A(ia,q) = 1;
        end
    end
    data(k).degree = degree(graph(A));
    data(k).A = A;
    [data(k).bins,data(k).binsize] = conncomp(graph(A));
    %Data organization for clm plotting
    %Array is [charge,z-location,l-per the cluster it belongs to, m-per the cluster it belongs to]
    cluster = zeros(length(A),5); %Pre-initializing array
    cluster(1:length(idxsa),1) = 1;
    cluster((length(idxsa)+1):end,1) = -1;
    cluster(1:length(idxsa),2) = xyz_cat(:,3);
    cluster((length(idxsa)+1):end,2) = data(k).com_a; %By our elegant preorganization of lammps file and our careful data
                                                      %analysis we have the following connection
    %Determining l-cations in the cluster and m-anions in the cluster
    kb=1; %Counting what bin we are looking at
    for lpm = data(k).binsize
        is = data(k).bins == kb; %indices of cations and anions involved in the cluster of binsize lpm = l+m
        l = sum(is(1:length(idxsa))); %Number of cations in the cluster
        s = length(unique([idxsw{is(1:length(idxsa))}])); %Number of solvent molecules bond to the cluster
        cluster(is,3:5) = cluster(is,3:5) + [l,lpm-l,s]; 
        kb = kb + 1;
    end
    data(k).cluster = cluster;
    data(k).typesofcluster = unique(strcat(num2str(cluster(:,3)),',',num2str(cluster(:,4)),',',num2str(cluster(:,5))),'rows');

    %Incerementing storage vector
    k=k+1;
    
end
fclose all;
end

function [r,dX,dY,dZ]=periodicDist(X,Y, L, zflag)

if nargin<3
    L=30;
    zflag=0;
end

if numel(L)==1
    L=repmat(L,1,3);
end

dX=bsxfun(@minus, X(:,1), Y(:,1)');
dY=bsxfun(@minus, X(:,2), Y(:,2)');
dZ=bsxfun(@minus, X(:,3), Y(:,3)');

dX=min(abs(dX), L(1)-abs(dX));
dY=min(abs(dY),L(2)-abs(dY));
if zflag
    dZ=min(abs(dZ),L(3)-abs(dZ));
end


r=sqrt(dX.^2+dY.^2+dZ.^2);
r=r';
end

function data = f_smoothdata_gaussian_int(data,ai,xvec,gw)

    nan_flag = 0;
    if any(isnan(data))%Here we have to cut more out, only case is in prob
        ai = ~isnan(data);
        nan_flag = 1;
    end
    data = data(ai);
    data = smoothdata(data,'gaussian',gw);

    if nan_flag == 0
        dataT = zeros(length(xvec),1); %Adding back in zeros
    else
        dataT = nan(length(xvec),1); %Adding back in nan
    end
    dataT(ai) = data;
    data = dataT;
end