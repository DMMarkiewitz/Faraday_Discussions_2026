function WiSE_EDL
%To get the average association constant we use the bulk values from the
%different surface charges we MD simulations for [0,0.01,0.05,0.1,0.2].
%Otherwise one only needs the 0.2 C/m^2 qs for the analysis in the paper 
clear all
close all
salt = 'LiTFSI';
qs = 0.2;%[0,0.01,0.05,0.1,0.2]%[0,0.01,0.05,0.1,0.15,0.2]%[0,0.01,0.05,0.1,0.15,0.2]; %[0,0.05,0.1,0.15,0.2]; % magnitude of surface charge density C/m^2
mm = 21; %molality
Lx = 32.624;
Ly = 32.624;
Lz = 265.248;
L = [Lx,Ly,Lz];
si = sqrt(1);
gw = 1;%if gw = 1 we are not smoothing the data
gw_p = 1;%if gw_p = 1 we are not smoothing the data
nbin = 650; %Number of bins
state = 1; %1 if we want to have the variations on 0 otherwise
sti = 1%0; %1 if we want it to be sticky otherwise 0
letterx = .6;
letters = 18;

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
if sti == 0
    fp = 5;
elseif sti == 1
    fp = 4;
end
fm = 3;
fs = 1;

%Ion association constants:
Lpm_p = 0;
Lpm_m = 0;
Lp0_p = 0;
Lp0_0 = 0;
L_tilde = 0;
%For stat's
a_Lpm_p = zeros(length(qs),1);
a_Lpm_m = zeros(length(qs),1);
a_Lp0_p = zeros(length(qs),1);
a_Lp0_0 = zeros(length(qs),1);
a_L_tilde = zeros(length(qs),1);
for i = 1 : length(qs)
    filePath = ['EDL_',num2str(qs(i)),'.xyz'];
    %[csalt,data]=readDumpFile(filePath,1,salt,mm,sep,cut_off,L);
    savename = ['S_EDL_paper26_sep_',num2str(sep),'_',num2str(qs(i)),'.mat'];
    %save(savename,'data','-v7.3')
    load(savename,'data')
    csalt = getconc(filePath,salt,mm,L);
    
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
        omeg = (xip+xim)*size(An_zpos,1)+size(Sol_zpos,1);
        Vbox = Lx*Ly*(Lz-2*Lx)*1e-30;%using potentially accessible region(Lzu-Lzl)
        rho_ao = nwo*size(Sol_zpos,1)/Vbox;
        rho_wo = nao*size(An_zpos,1)/Vbox;
        v0 = Vbox/omeg; %m^3 volume backcalculated from vbox and Omega
        T = 300; 
        Na = 6.022e23; 
        kB = 1.38e-23;  
        e = 1.6e-19;
        beta = 1/(kB*T);
        eps_r = 10.1; %So our debye length used to non-dimensionalize is same as our theory
        eps=eps_r*8.85e-12;
        phi_0_b = size(Sol_zpos,1)/omeg;
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
        mean_Li_c_an = mean_Li_c_an + mean(data(n).coord_li_an(ceil(end/3):ceil(2*end/3),2))/n_length;
        mean_Li_c_w = mean_Li_c_w + mean(data(n).coord_li_w(ceil(end/3):ceil(2*end/3),2))/n_length;
    end

    %Alternative smoothing
    Li_coord_ans = f_smoothdata_gaussian_int(Li_coord_an,active_index,x,gw);
    Li_coord_ws = f_smoothdata_gaussian_int(Li_coord_w,active_index,x,gw);

    %Calculations for association probabilities  
    %As there is values associated with the positions we need to first
    %parse our data and do our summing ourself
    pO = zeros(length(x),1); %Preinitialize pO-array, coords are the x bins coords
    Op = pO; %Preinitialize Op-array, coords are the x bins coords
    pm = pO; %Preinitialize pm-array, coords are the x bins coords
    mp = pO; %Preinitialize mp-array, coords are the x bins coords
    nc = pO; %Preinitialize nc-array, coords are the x bins coords
    na = pO; %Preinitialize na-array, coords are the x bins coords
    ns = pO; %Preinitialize ns-array, coords are the x bins coords
    n_length = length(data); %How many time steps we are averaging our data over
    n_cat = pO; %Tracks how many timesteps the box is truely active for i.e. have one particle of their class present in
    n_an = pO; %Tracks how many timesteps the box is truely active for i.e. have one particle of their class present in
    n_sol = pO; %Tracks how many timesteps the box is truely active for i.e. have one particle of their class present in
    for n = 1:n_length
        Ncat = histcounts(data(n).com_c,ed);
        Nan = histcounts(data(n).com_a,ed);
        Nsol = histcounts(data(n).xyz_ow(:,3),ed);
        nc = nc + Ncat'/n_length;
        na = na + Nan'/n_length;
        ns = ns + Nsol'/n_length;
        for ii = 1:length(x) % for loop running over the bins
            inx1 = (ed(ii) <= data(n).pOassoc(:,1)) & (data(n).pOassoc(:,1) <= ed(ii+1));
            inx2 = (ed(ii) <= data(n).Opassoc(:,1)) & (data(n).Opassoc(:,1) <= ed(ii+1));
            inx3 = (ed(ii) <= data(n).pmassoc(:,1)) & (data(n).pmassoc(:,1) <= ed(ii+1));
            inx4 = (ed(ii) <= data(n).mpassoc(:,1)) & (data(n).mpassoc(:,1) <= ed(ii+1));
            
            if Ncat(ii) ~= 0
                pO(ii) = pO(ii) + sum(data(n).pOassoc(inx1,2))/Ncat(ii); %Summing association number in this bin and weigthing it unifromly amond the timesteps and dividing by fp and the number of the molecules in this region
                pm(ii) = pm(ii) + sum(data(n).pmassoc(inx3,2))/Ncat(ii); %Summing association number in this bin and weigthing it unifromly amond the timesteps and dividing by fp and the number of the molecules in this region
                n_cat(ii) = n_cat(ii) + 1;
            end

            if Nan(ii) ~= 0
                mp(ii) = mp(ii) + sum(data(n).mpassoc(inx4,2))/Nan(ii); %Summing association number in this bin and weigthing it unifromly amond the timesteps and dividing by fm and the number of the molecules in this region
                n_an(ii) = n_an(ii) + 1;
            end
                 
            if Nsol(ii) ~= 0
                Op(ii) = Op(ii) + sum(data(n).Opassoc(inx2,2))/Nsol(ii); %Summing association number in this bin and weigthing it unifromly amond the timesteps and dividing by fs and the number of the molecules in this region
                n_sol(ii) = n_sol(ii) + 1;
            end

        end
    end

    %Calculating the pij(x) temporally averaged
    TpO = pO./(fp*n_cat);
    Tpm = pm./(fp*n_cat);
    mp = mp./(fm*n_an);
    Op = Op./(fs*n_sol);

    if sti ~= 1
        %For Non-sticky values when fp = 5 and non-renormalization occurs
        pO = TpO;
        pm = Tpm;
    else
        %For the sticky case
        pO = TpO*fp./(mean_Li_c_an + mean_Li_c_w);
        pm = Tpm*fp./(mean_Li_c_an + mean_Li_c_w);
    end

    %Getting our variance for pij and ni
    vpO = zeros(length(x),1); %Preinitialize pO-array, coords are the x bins coords
    vOp = vpO; %Preinitialize Op-array, coords are the x bins coords
    vpm = vpO; %Preinitialize pm-array, coords are the x bins coords
    vmp = vpO; %Preinitialize mp-array, coords are the x bins coords
    vnc = vpO; %Preinitialize nc-array, coords are the x bins coords
    vna = vpO; %Preinitialize na-array, coords are the x bins coords
    vns = vpO; %Preinitialize ns-array, coords are the x bins coords
    n_length = length(data); %How many time steps we are averaging our data over
    n_cat = vpO; %Tracks how many timesteps the box is truely active for i.e. have one particle of their class present in
    n_an = vpO; %Tracks how many timesteps the box is truely active for i.e. have one particle of their class present in
    n_sol = vpO; %Tracks how many timesteps the box is truely active for i.e. have one particle of their class present in
    for n = 1:n_length
        Ncat = histcounts(data(n).com_c,ed);
        Nan = histcounts(data(n).com_a,ed);
        Nsol = histcounts(data(n).xyz_ow(:,3),ed);
        vnc = vnc + (Ncat'-nc).^2/(n_length-1);
        vna = vna + (Nan'-na).^2/(n_length-1);
        vns = vns + (Nsol'-ns).^2/(n_length-1);
        for ii = 1:length(x) % for loop running over the bins
            inx1 = (ed(ii) <= data(n).pOassoc(:,1)) & (data(n).pOassoc(:,1) <= ed(ii+1));
            inx2 = (ed(ii) <= data(n).Opassoc(:,1)) & (data(n).Opassoc(:,1) <= ed(ii+1));
            inx3 = (ed(ii) <= data(n).pmassoc(:,1)) & (data(n).pmassoc(:,1) <= ed(ii+1));
            inx4 = (ed(ii) <= data(n).mpassoc(:,1)) & (data(n).mpassoc(:,1) <= ed(ii+1));
            if Ncat(ii) ~= 0
                if sti ~= 1
                    % For Non sticky case
                    vpO(ii) = vpO(ii) + (sum(data(n).pOassoc(inx1,2))/(fp*Ncat(ii))-pO(ii))^2; %Summing association number in this bin and weigthing it unifromly amond the timesteps and dividing by fp and the number of the molecules in this region
                    vpm(ii) = vpm(ii) + (sum(data(n).pmassoc(inx3,2))/(fp*Ncat(ii))-pm(ii))^2; %Summing association number in this bin and weigthing it unifromly amond the timesteps and dividing by fp and the number of the molecules in this region
                else
                    % For Sticky Case
                    vpO(ii) = vpO(ii) + (sum(data(n).pOassoc(inx1,2))/((mean_Li_c_an + mean_Li_c_w)*Ncat(ii))-pO(ii))^2;
                    vpm(ii) = vpm(ii) + (sum(data(n).pmassoc(inx3,2))/((mean_Li_c_an + mean_Li_c_w)*Ncat(ii))-pm(ii))^2;
                end
                n_cat(ii) = n_cat(ii) + 1;
            end

            if Nan(ii) ~= 0
                vmp(ii) = vmp(ii) + (sum(data(n).mpassoc(inx4,2))/(fm*Nan(ii))-mp(ii))^2; %Summing association number in this bin and weigthing it unifromly amond the timesteps and dividing by fm and the number of the molecules in this region
                n_an(ii) = n_an(ii) + 1;
            end
                 
            if Nsol(ii) ~= 0
                vOp(ii) = vOp(ii) + (sum(data(n).Opassoc(inx2,2))/Nsol(ii)-Op(ii))^2; %Summing association number in this bin and weigthing it unifromly amond the timesteps and dividing by fs and the number of the molecules in this region
                n_sol(ii) = n_sol(ii) + 1;
            end

        end
    end

    %Calculating the sd of ppij(x)
    sdpO = sqrt(vpO./max(0,n_cat-1));
    sdpm = sqrt(vpm./max(0,n_cat-1));
    sdmp = sqrt(vmp./max(0,n_an-1));
    sdOp = sqrt(vOp./max(0,n_sol-1));

    %Calculating the sd of Ni(x)
    sdnc = sqrt(vnc);
    sdna = sqrt(vna);
    sdns = sqrt(vns);

    %Calculating the other sigma with propogation of error
    sdsxin = sqrt((xip*sdnc).^2+(xim*sdna).^2+sdns.^2);
    sdgel = pm.*mp.*sqrt((sdpm./pm).^2+(sdmp./mp).^2);
    sdphic = nc.*xip./(nc*xip+na*xim+ns).*sqrt((sdnc./nc).^2+(sdsxin./(nc*xip+na*xim+ns)).^2);
    sdphia = na.*xim./(nc*xip+na*xim+ns).*sqrt((sdna./na).^2+(sdsxin./(nc*xip+na*xim+ns)).^2);
    sdphis = ns./(nc*xip+na*xim+ns).*sqrt((sdns./ns).^2+(sdsxin./(nc*xip+na*xim+ns)).^2);

    %Correcting for cases where the mean value is zero
    sdsxin(isnan(sdsxin)) = 0;
    sdgel(isnan(sdgel)) = 0;
    sdphic(isnan(sdphic)) = 0;
    sdphia(isnan(sdphia)) = 0;
    sdphis(isnan(sdphis)) = 0;

    %Smoothed Versions
    %Generating Smoothed versions
    pms = f_smoothdata_gaussian_int(pm,active_index,x,gw_p);
    mps = f_smoothdata_gaussian_int(mp,active_index,x,gw_p);
    pOs = f_smoothdata_gaussian_int(pO,active_index,x,gw_p);
    Ops = f_smoothdata_gaussian_int(Op,active_index,x,gw_p);

    %New distibution calculations
    %Finding all the clusters we see in all timesteps
    clusterlist = [0,0,1]; 
    for ia = 1:n
       clusterlist = cat(1,clusterlist,str2double(split(string(data(ia).typesofcluster),',')));
    end
    tot_typesofcluster = unique(clusterlist,'rows'); %Columns are l,m,s

    %Preinitializing cluster info storage
    Nlms = zeros(size(tot_typesofcluster,1),length(x)); %Stores the clm info rows corespond to order in types of cluster and columns are varying z with positions shown in z

    %Obtaining time averaged distribution Clm(z)
    for di = 1:length(data)
        %Obtaining current distribution Clm(z) at fixed t
        for cti = 1:size(tot_typesofcluster,1)
            l = tot_typesofcluster(cti,1);
            m = tot_typesofcluster(cti,2);
            s = tot_typesofcluster(cti,3);
            if cti~=1
                ap = all(data(di).cluster(:,3:5) == [l,m,s],2); %Active particle indices
                [TNlms,~] = histcounts(data(di).cluster(ap,2),ed); %Total number of lm cations in delta z
                Nlms(cti,:) = Nlms(cti,:) + TNlms/((l+m)*length(data));
            else
                [TNlms,~] = histcounts(data(di).com_s(~data(di).label_s),ed);
                Nlms(cti,:) = Nlms(cti,:) + TNlms/(length(data));
            end
        end
    end

    %Turning Nlms to clms
    clms = Nlms/omeg;
    
    %Pre-initializing
    phi_lms = Nlms;
    
    %Getting all phi_lms        
    for cti = 1:length(tot_typesofcluster)
        Nlms(cti,:) = clms(cti,:)*(tot_typesofcluster(cti,:)*[xip;xim;1]);
    end

    %Extracting Nlms's of interest if they exist
    iN_100 = zeros(1,length(x));
    iN_010 = zeros(1,length(x));
    iN_001 = zeros(1,length(x));
    iN_10x = zeros(1,length(x));
    iN_11x = zeros(1,length(x));
    iN_g11x = zeros(1,length(x));

    for cti = 1:length(tot_typesofcluster)
        if all(tot_typesofcluster(cti,:)==[1,0,0])
            iN_100 = Nlms(cti,:);
        elseif all(tot_typesofcluster(cti,:)==[0,1,0])
            iN_010 = Nlms(cti,:);
        elseif all(tot_typesofcluster(cti,:)==[0,0,1])
            iN_001 = Nlms(cti,:);
        elseif all(tot_typesofcluster(cti,1:2)==[1,0])
            iN_10x = iN_10x + Nlms(cti,:);
        elseif all(tot_typesofcluster(cti,1:2)==[1,1])
            iN_11x = iN_11x + Nlms(cti,:);
        else
            iN_g11x = iN_g11x + Nlms(cti,:);
        end
    end
    
    %Need to remormalize
    N_norm = 1./(iN_100+iN_010+iN_001+iN_10x+iN_11x+iN_g11x);

    x_10x = (iN_10x+iN_100).*N_norm;
    x_010 = iN_010.*N_norm;
    x_001 = iN_001.*N_norm;
    x_Agg = (iN_11x+iN_g11x).*N_norm;

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

    %Finding the region of expected model breakdown
    [gl,~] = max([min(x(nc~=0)),min(x(na~=0)),min(x(ns~=0))]);
    [gu,~] = min([max(x(nc~=0)),max(x(na~=0)),max(x(ns~=0))]);
    xminvf = xmin;
    xmaxvf = xmax;

    %Plotting Coordination
    figure
    plot((x(act_i)-xmin)*AA2kappa,Li_coord_ans(act_i)*mean_Li_c_an/mean(Li_coord_ans(ceil(end/3):ceil(2*end/3))),'b-',(x(act_i)-xmin)*AA2kappa,Li_coord_ws(act_i)*mean_Li_c_w/mean(Li_coord_ws(ceil(end/3):ceil(2*end/3))),'k--',(x(act_i)-xmin)*AA2kappa,Li_coord_ws(act_i)*mean_Li_c_w/mean(Li_coord_ws(ceil(end/3):ceil(2*end/3)))+Li_coord_ans(act_i)*mean_Li_c_an/mean(Li_coord_ans(ceil(end/3):ceil(2*end/3))),'r:','linewidth',1.5)
    ylim([0,8])
    xlim([0,(max(x(act_i)-xmin)*AA2kappa)])
    legend({['Li$^+$ (Anion), Bulk Avg.\ = ',num2str(round(mean_Li_c_an,1),'%.1f')],['Li$^+$ (H$_2$O), Bulk Avg.\ = ',num2str(round(mean_Li_c_w,1),'%.1f')],['Li$^+$ (Total), Bulk Avg.\ = ',num2str(round(mean_Li_c_an,1)+round(mean_Li_c_w,1),'%.1f')]},'interpreter','latex', 'location','north')
    xlabel('$\kappa x$','interpreter','latex')
    ylabel('Li$^+$ Coordination Number','interpreter','latex')
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 16);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'color', 'none');
    set(g,'InvertHardcopy','on');
    saveas(gcf,strcat('Cordination_Number_For_',num2str(mm),'m_q',strrep(num2str(qs(i)),'.','_'),'_bin_',num2str(nbin),'_sti_',num2str(sti)),'epsc')
    saveas(gcf,strcat('Cordination_Number_For_',num2str(mm),'m_q',strrep(num2str(qs(i)),'.','_'),'_bin_',num2str(nbin),'_sti_',num2str(sti)),'jpeg')
    savefig(strcat('Cordination_Number_For_',num2str(mm),'m_q',strrep(num2str(qs(i)),'.','_'),'_bin_',num2str(nbin),'_sti_',num2str(sti),'.fig'))

    %Making tiled plots
    figure('Renderer', 'painters', 'units','inches','Position',[.01 .01 9 9])
    t = tiledlayout(4,2,'TileSpacing','tight','padding','compact')

    %Plotting Volume Fractions
    %Positive side
    nexttile
    hold on
    plot((x(act_i)-xmin)*AA2kappa,nc(act_i)'*xip./(nc(act_i)'*xip+na(act_i)'*xim+ns(act_i)'),'r-','linewidth',1.5)
    plot((x(act_i)-xmin)*AA2kappa,na(act_i)'*xim./(nc(act_i)'*xip+na(act_i)'*xim+ns(act_i)'),'b--','linewidth',1.5)
    plot((x(act_i)-xmin)*AA2kappa,ns(act_i)'./(nc(act_i)'*xip+na(act_i)'*xim+ns(act_i)'),'k:','linewidth',1.5)
    xregion(0,(gl-xminvf)*AA2kappa)
    xlim([0,30])
    ylim([-0.1,1.29])
    box on
    leg=legend({'Li$^+$','TFSI$^-$','H$_2$O'},'location','north','interpreter','latex','Orientation','horizontal','box','off');
    leg.ItemTokenSize = [18,18];
    ylabel('$\bar{\phi}_i$','interpreter','latex')
    title('')
    text(letterx,1.15,'\bf{a)}','interpreter','latex','FontSize',letters)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    ax = gca;
    ax.FontSize = 16;
    ax.XTick = 0:5:30;
    ax.YTick = 0:.2:1;
    set(gca,'XTickLabel',[]);
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 16);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'color', 'none');
    set(g,'InvertHardcopy','on');
    
    %Negative side
    nexttile
    hold on
    plot([0,(x(act_i)-xmin)*AA2kappa],[0,fliplr(nc(act_i)')*xip./(fliplr(nc(act_i)')*xip+fliplr(na(act_i)')*xim+fliplr(ns(act_i)'))],'r-','linewidth',1.5)
    plot([0,(x(act_i)-xmin)*AA2kappa],[0,fliplr(na(act_i)')*xim./(fliplr(nc(act_i)')*xip+fliplr(na(act_i)')*xim+fliplr(ns(act_i)'))],'b--','linewidth',1.5)
    plot([0,(x(act_i)-xmin)*AA2kappa],[0,fliplr(ns(act_i)')./(fliplr(nc(act_i)')*xip+fliplr(na(act_i)')*xim+fliplr(ns(act_i)'))],'k:','linewidth',1.5)
    xregion(0,(xmaxvf-gu)*AA2kappa)
    xlim([0,30])
    ylim([-0.1,1.29])
    box on
    leg=legend({'Li$^+$','TFSI$^-$','H$_2$O'},'location','north','interpreter','latex','Orientation','horizontal','box','off');
    leg.ItemTokenSize = [18,18];
    ylabel(' ')
    text(letterx,1.15,'\bf{e)}','interpreter','latex','FontSize',letters)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    ax = gca;
    ax.FontSize = 16;
    ax.XTick = 0:5:30;
    ax.YTick = 0:.2:1;
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 16);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'color', 'none');
    set(g,'InvertHardcopy','on');

    %Plotting Cluster Volume Fractions
    %Positive side
    nexttile
    hold on
    plot((x(act_i)-xmin)*AA2kappa,x_10x(act_i),'r-','linewidth',1.5)
    plot((x(act_i)-xmin)*AA2kappa,x_010(act_i),'b--','linewidth',1.5)
    plot((x(act_i)-xmin)*AA2kappa,x_001(act_i),'k:','linewidth',1.5)
    plot((x(act_i)-xmin)*AA2kappa,x_Agg(act_i),'m-.','linewidth',1.5)
    xregion(0,(gl-xminvf)*AA2kappa)
    xlim([0,30])
    ylim([-0.1,1.29])
    box on
    leg=legend({'$\bar{\phi}_{10x}$','$\bar{\phi}_{010}$','$\bar{\phi}_{001}$','$\bar{\phi}_{Agg}$'},'location','north','interpreter','latex','Orientation','horizontal','box','off');
    leg.ItemTokenSize = [18,18];
    ylabel('$\bar{\phi}_{lms}$','interpreter','latex')
    text(letterx,1.15,'\bf{b)}','interpreter','latex','FontSize',letters)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    ax = gca;
    ax.FontSize = 16;
    ax.XTick = 0:5:30;
    ax.YTick = 0:.2:1;
    set(gca,'XTickLabel',[]);
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 16);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'color', 'none');
    set(g,'InvertHardcopy','on');
    
    %Negative side
    nexttile
    hold on   
    plot([0,(x(act_i)-xmin)*AA2kappa],[0,fliplr(x_10x(act_i))],'r-','linewidth',1.5)
    plot([0,(x(act_i)-xmin)*AA2kappa],[0,fliplr(x_010(act_i))],'b--','linewidth',1.5)
    plot([0,(x(act_i)-xmin)*AA2kappa],[0,fliplr(x_001(act_i))],'k:','linewidth',1.5)
    plot([0,(x(act_i)-xmin)*AA2kappa],[0,fliplr(x_Agg(act_i))],'m-.','linewidth',1.5)
    xregion(0,(xmaxvf-gu)*AA2kappa)
    xlim([0,30])
    ylim([-0.1,1.29])
    box on
    leg=legend({'$\bar{\phi}_{10x}$','$\bar{\phi}_{010}$','$\bar{\phi}_{001}$','$\bar{\phi}_{Agg}$'},'location','north','interpreter','latex','Orientation','horizontal','box','off');
    leg.ItemTokenSize = [18,18];
    text(letterx,1.15,'\bf{f)}','interpreter','latex','FontSize',letters)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    ax = gca;
    ax.FontSize = 16;
    ax.XTick = 0:5:30;
    ax.YTick = 0:.2:1;
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 16);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'color', 'none');
    set(g,'InvertHardcopy','on');

    %Plotting Association Probabilities
    %Positive side
    nexttile
    hold on
    plot((x(act_i)-xmin)*AA2kappa,pms(act_i),'r-',(x(act_i)-xmin)*AA2kappa,pOs(act_i),'r--','linewidth',1.5)
    plot((x(act_i)-xmin)*AA2kappa,mps(act_i),'b:',(x(act_i)-xmin)*AA2kappa,Ops(act_i),'k-.','linewidth',1.5)
    xregion(0,(gl-xminvf)*AA2kappa)
    xlim([0,30])
    ylim([-0.1,1.39])
    box on
    leg=legend({'$\bar{p}_{+-}$','$\bar{p}_{+0}$','$\bar{p}_{-+}$','$\bar{p}_{0+}$'},'location','north','interpreter','latex','Orientation','horizontal','box','off');
    leg.ItemTokenSize = [18,18];
    ylabel('$\bar{p}_{ij}$','interpreter','latex')
    text(letterx,1.25,'\bf{c)}','interpreter','latex','FontSize',letters)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    ax = gca;
    ax.FontSize = 16;
    ax.XTick = 0:5:30;
    ax.YTick = 0:.2:1.2;
    set(gca,'XTickLabel',[]);
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 16);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'color', 'none');
    set(g,'InvertHardcopy','on');

    %Negative side
    nexttile
    hold on
    plot((x(act_i)-xmin)*AA2kappa,fliplr(pms(act_i)'),'r-',(x(act_i)-xmin)*AA2kappa,fliplr(pOs(act_i)'),'r--','linewidth',1.5)
    plot((x(act_i)-xmin)*AA2kappa,fliplr(mps(act_i)'),'b:',(x(act_i)-xmin)*AA2kappa,fliplr(Ops(act_i)'),'k-.','linewidth',1.5)
    xregion(0,(xmaxvf-gu)*AA2kappa)
    xlim([0,30])
    ylim([-0.1,1.39])
    box on
    leg=legend({'$\bar{p}_{+-}$','$\bar{p}_{+0}$','$\bar{p}_{-+}$','$\bar{p}_{0+}$'},'location','north','interpreter','latex','Orientation','horizontal','box','off');
    leg.ItemTokenSize = [18,18];
    text(letterx,1.25,'\bf{g)}','interpreter','latex','FontSize',letters)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    ax = gca;
    ax.FontSize = 16;
    ax.XTick = 0:5:30;
    ax.YTick = 0:.2:1.2;
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 16);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'color', 'none');
    set(g,'InvertHardcopy','on');

    %Plotting Proximity to gelation
    %Positive side
    nexttile
    hold on
    plot((x(act_i)-xmin)*AA2kappa,1/((fp-1)*(fm-1))*ones(length(x(act_i)),1),'k--',(x(act_i)-xmin)*AA2kappa,pms(act_i).*mps(act_i),'r-','linewidth',1.5)
    shplot(state,(x(act_i)-xmin)*AA2kappa,pms(act_i).*mps(act_i),sdgel(act_i),'r-')
    plot((x(act_i)-xmin)*AA2kappa,1/((fp-1)*(fm-1))*ones(length(x(act_i)),1),'k--','linewidth',1.5)
    xregion(0,(gl-xminvf)*AA2kappa)
    xlim([0,30])
    ylim([-0.05,.675])
    box on
    leg=legend('$((f_+-1)(f_--1))^{-1}$','location','north','interpreter','latex','Orientation','horizontal','box','off');
    leg.ItemTokenSize = [18,18];
    xlabel('$\kappa x$','interpreter','latex')
    ylabel('$\bar{p}_{+-}\bar{p}_{-+}$','interpreter','latex')
    text(letterx,0.6,'\bf{d)}','interpreter','latex','FontSize',letters)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    ax = gca;
    ax.FontSize = 16;
    ax.XTick = 0:5:30;
    ax.YTick = 0:.1:.8;
    g=gcf;
    g.Renderer='painters';
    set(gca, 'Xcolor', 'k');
    set(gca, 'Ycolor', 'k');
    set(gca, 'FontSize', 16);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'Layer', 'Top');
    set(gca, 'color', 'none');
    set(g,'InvertHardcopy','on');

    %Negative side
    nexttile
    hold on
    plot((x(act_i)-xmin)*AA2kappa,1/((fp-1)*(fm-1))*ones(length(x(act_i)),1),'k--',(x(act_i)-xmin)*AA2kappa,fliplr((pms(act_i).*mps(act_i))'),'r-','linewidth',1.5)
    shplot(state,(x(act_i)-xmin)*AA2kappa,fliplr((pms(act_i).*mps(act_i))'),fliplr(sdgel(act_i)'),'r-')
    plot((x(act_i)-xmin)*AA2kappa,1/((fp-1)*(fm-1))*ones(length(x(act_i)),1),'k--','linewidth',1.5)
    xregion(0,(xmaxvf-gu)*AA2kappa)
    xlim([0,30])        
    ylim([-0.05,.675])
    box on
    leg=legend('$((f_+-1)(f_--1))^{-1}$','location','north','interpreter','latex','Orientation','horizontal','box','off');
    leg.ItemTokenSize = [18,18];
    xlabel('$\kappa x$','interpreter','latex')
    text(letterx,0.6,'\bf{h)}','interpreter','latex','FontSize',letters)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'DefaultTextInterpreter', 'latex')
    ax = gca;
    ax.FontSize = 16;
    ax.XTick = 0:5:30;
    ax.YTick = 0:.1:.8;
    set(gca,'YTickLabel',[]);
    set(gca, 'Layer', 'top');
    set(gca, 'LineWidth', 1.5);
    set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','on');
    
    %Saving the panel figure
    exportgraphics(t,strcat(strrep(strcat('Plots_For_',num2str(mm),'m_q',num2str(qs(i)),'_bin_',num2str(nbin),'_sti_',num2str(sti),'_er_',num2str(state)),'.','_'),'.eps'),'Resolution',600)
    exportgraphics(t,strcat(strrep(strcat('Plots_For_',num2str(mm),'m_q',num2str(qs(i)),'_bin_',num2str(nbin),'_sti_',num2str(sti),'_er_',num2str(state)),'.','_'),'.jpeg'),'Resolution',600)
    savefig(strcat(strrep(strcat('Plots_For_',num2str(mm),'m_q',num2str(qs(i)),'_bin_',num2str(nbin),'_sti_',num2str(sti),'_er_',num2str(state)),'.','_'),'.fig'))
    
    %Calculating bulk probabilities
    pm_b = pm(length(Lzl:-(Lzu-Lzl)/(nbin-1):0):length(unique([flip(Lzl:-(Lzu-Lzl)/(nbin-1):0) linspace(Lzl,Lzu,nbin)])));
    pm_b = mean(pm_b(ceil(nbin/3):floor(2*nbin/3)));
    mp_b = mp(length(Lzl:-(Lzu-Lzl)/(nbin-1):0):length(unique([flip(Lzl:-(Lzu-Lzl)/(nbin-1):0) linspace(Lzl,Lzu,nbin)])));
    mp_b = mean(mp_b(ceil(nbin/3):floor(2*nbin/3)));
    pO_b = pO(length(Lzl:-(Lzu-Lzl)/(nbin-1):0):length(unique([flip(Lzl:-(Lzu-Lzl)/(nbin-1):0) linspace(Lzl,Lzu,nbin)])));
    pO_b = mean(pO_b(ceil(nbin/3):floor(2*nbin/3)));
    Op_b = Op(length(Lzl:-(Lzu-Lzl)/(nbin-1):0):length(unique([flip(Lzl:-(Lzu-Lzl)/(nbin-1):0) linspace(Lzl,Lzu,nbin)])));
    Op_b = mean(Op_b(ceil(nbin/3):floor(2*nbin/3)));

    %Calculating bulk phi's
    phi_p = nc*xip./(nc*xip+na*xim+ns);
    phi_m = na*xim./(nc*xip+na*xim+ns);
    phi_0 = ns./(nc*xip+na*xim+ns);

    phi_p_b = phi_p(length(Lzl:-(Lzu-Lzl)/(nbin-1):0):length(unique([flip(Lzl:-(Lzu-Lzl)/(nbin-1):0) linspace(Lzl,Lzu,nbin)])));
    phi_p_b = mean(phi_p_b(ceil(nbin/3):floor(2*nbin/3)));
    phi_m_b = phi_m(length(Lzl:-(Lzu-Lzl)/(nbin-1):0):length(unique([flip(Lzl:-(Lzu-Lzl)/(nbin-1):0) linspace(Lzl,Lzu,nbin)])));
    phi_m_b = mean(phi_m_b(ceil(nbin/3):floor(2*nbin/3)));
    phi_0_b = phi_0(length(Lzl:-(Lzu-Lzl)/(nbin-1):0):length(unique([flip(Lzl:-(Lzu-Lzl)/(nbin-1):0) linspace(Lzl,Lzu,nbin)])));
    phi_0_b = mean(phi_0_b(ceil(nbin/3):floor(2*nbin/3)));

    if sti ~= 1
        %For Non-sticky values when fp = 5 and non-renormalization occurs
        a_Lpm_p(i) = mp_b*xip/(fp*phi_p_b*(1-pm_b-pO_b)*(1-mp_b));
        a_Lpm_m(i) = pm_b*xim/(fm*phi_m_b*(1-pm_b-pO_b)*(1-mp_b));
        a_Lp0_p(i) = Op_b*xip/(fp*phi_p_b*(1-pm_b-pO_b)*(1-Op_b));
        a_Lp0_0(i) = pO_b/(phi_0_b*(1-pm_b-pO_b)*(1-Op_b));
    
        Lpm_p = Lpm_p + a_Lpm_p(i)/length(qs);
        Lpm_m = Lpm_m + a_Lpm_m(i)/length(qs);
        Lp0_p = Lp0_p + a_Lp0_p(i)/length(qs);
        Lp0_0 = Lp0_0 + a_Lp0_0(i)/length(qs);
    else
        a_L_tilde(i) = mp_b*(1-Op_b)/(Op_b*(1-mp_b));
    
        L_tilde = L_tilde + a_L_tilde(i)/length(qs);
    end
end

%Reporting the current values predictions
if sti ~= 1
    %Non-sticky case
    Lpm_p
    Lpm_m
    Lpm = (Lpm_p + Lpm_m)/2
    s_Lpm_p = sqrt(sum((a_Lpm_p-Lpm).^2)/(length(qs)-1))
    s_Lpm_m = sqrt(sum((a_Lpm_m-Lpm).^2)/(length(qs)-1))
    Lp0_p
    Lp0_0 
    Lp0 = (Lp0_p + Lp0_0)/2
    s_Lp0_p = sqrt(sum((a_Lp0_p-Lp0).^2)/(length(qs)-1))
    s_Lp0_0 = sqrt(sum((a_Lp0_0-Lp0).^2)/(length(qs)-1))
    L_tilde = Lpm/Lp0
else
    %Sticky case
    L_tilde
    s_L_tilde = sqrt(sum((a_L_tilde-L_tilde).^2)/(length(qs)-1))
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

function varargout = shplot(state,x,y,e,varargin)
%% H = shplot(x,y,e,opt,varargin)
% This function creates a shaded error-bar plot. It is a lightweight
% implementation of 'shadedErrorBar' by Rob Campbell intended
% to be used as a local function in other scripts.
%
% H is a structure of handles to parts of the plot (line path upper lower).
%
% EXAMPLE
%  n = 18; x = linspace(-1,1,n^2); y = peaks(n); y = y(:);
%  for j = 1:100, Y(:,j) = y + randn(n^2,1); end;
%  opt = {'Color', [1 .39 .3], 'Marker', '.'};
%  H = shplot(x,mean(Y,2),6*std(Y,[],2)./sqrt(size(Y,2)), opt{:});
%  grid on; box on; xlim([-1 1]);
%  legend([H.line H.patch],'<Y>','6\sigma_{m}(Y)');
%
% Written by Marcin Konowalczyk
% Timmel Group @ Oxford University
% Edited by Daniel M. Markiewitz @ MIT
x = x(:); y = y(:); e = abs(e(:));
isheld = ishold; if ~isheld; cla; hold on; end
H.line = plot(x,y,varargin{:},'linewidth',1.5);
eu = min(1,y + e); el = max(0,y - e);
col = 0.15*get(H.line,'color') + 0.85; 
ecol = 3*col-2;
lst = get(H.line,'LineStyle');
set(gcf,'renderer','painters');
if state==1
    a=0.6;
else
    a=0;
end
H.patch = patch([x(~isnan(y)) ;flipud(x(~isnan(y)))],[el(~isnan(y)) ;flipud(eu(~isnan(y)))],1,'facecolor',col,'edgecolor','none','facealpha',a);
uistack(H.line,'top');
if ~isheld; hold off; end
if nargout > 0; varargout = {H}; end
end