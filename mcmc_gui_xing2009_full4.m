function mcmc_gui_xing2009_full
% ============================================================
% 모델링 회의체 화이팅!
% 모델링 회의체 화이팅!!
% FULL Xing 2009 6-ODE MCMC GUI (Improved Version)
% - Excel import (header auto-detect + unit convert)
% - Full 6-variable ODE (Xing 2009 기반식)
% - Stoichiometric 6개 자동 추정
% - MCMC 12개 파라미터 + Adaptive Proposal auto-tune
% - Likelihood weighting 개선 (변수별 σ 자동 스케일링)
% - Batch별 실험값 + Batch별 예측곡선 + 평균선 출력
% ============================================================

clc; close all;

%% ---------------------- GUI ----------------------
f = figure('Name','Full Xing2009 6-ODE MCMC','Color','w',...
    'Position',[200 80 850 700]);

uicontrol(f,'Style','text','String','엑셀 파일 (각 시트=1 batch):',...
    'BackgroundColor','w','HorizontalAlignment','left',...
    'Position',[20 650 400 25]);

hFile = uicontrol(f,'Style','edit','BackgroundColor','w',...
    'Position',[20 620 500 28],'String','');

uicontrol(f,'Style','pushbutton','String','Load Excel',...
    'Position',[540 620 130 28],'BackgroundColor',[0.2 0.6 1],...
    'ForegroundColor','w','Callback',@load_excel);

uicontrol(f,'Style','text','String','MCMC Iteration:',...
    'BackgroundColor','w','HorizontalAlignment','left',...
    'Position',[20 570 200 25]);

hIter = uicontrol(f,'Style','edit','BackgroundColor','w',...
    'Position',[20 545 130 28],'String','2000');

hRun = uicontrol(f,'Style','pushbutton','String','Run MCMC',...
    'Position',[260 500 190 45],'BackgroundColor',[0.2 0.6 1],...
    'ForegroundColor','w','FontWeight','bold','Callback',@run_mcmc);

hStatus = uicontrol(f,'Style','text','String','엑셀을 먼저 불러오세요.',...
    'Position',[20 470 800 28],'HorizontalAlignment','center',...
    'FontWeight','bold','BackgroundColor','w');

batches = {};   % loaded Excel batches

%% ============================================================
% 1) LOAD EXCEL
% ============================================================
function load_excel(~,~)

    [file,path] = uigetfile({'*.xlsx;*.xls'});
    if isequal(file,0), return; end

    filepath = fullfile(path,file);
    set(hFile,'String',filepath);

    [~, sheets] = xlsfinfo(filepath);
    N = numel(sheets);

    B = {};
    for s = 1:N
        T = readtable(filepath,'Sheet',sheets{s},...
            'PreserveVariableNames',true);

        headers = T.Properties.VariableNames;

        dat = struct('t',[],'x',[],'via',[],...
            'glc',[],'gln',[],'lac',[],'nh4',[]);

        units = struct('t','','x','','via','',...
            'glc','','gln','','lac','','nh4','');

        for c = 1:numel(headers)
            raw = headers{c};
            h = lower(regexprep(raw,'\s+',''));
            tok = regexp(raw,'\((.*?)\)','tokens','once');
            u = '';
            if ~isempty(tok), u = tok{1}; end

            vec = toNum(T{:,c});

            if contains(h,'time')||contains(h,'hour')
                dat.t = vec; units.t = u;
            elseif contains(h,'vcd')||contains(h,'cell')
                dat.x = vec; units.x = u;
            elseif contains(h,'via')||contains(h,'viab')
                dat.via = vec; units.via = u;
            elseif contains(h,'glc')||contains(h,'glucose')
                dat.glc = vec; units.glc = u;
            elseif contains(h,'gln')||contains(h,'glutamine')
                dat.gln = vec; units.gln = u;
            elseif contains(h,'lac')||contains(h,'lactate')
                dat.lac = vec; units.lac = u;
            elseif contains(h,'nh4')||contains(h,'amm')
                dat.nh4 = vec; units.nh4 = u;
            end
        end

        if isempty(dat.t), continue; end

        dat = unitConvert(dat,units);

        mask = ~isnan(dat.t);
        t = dat.t(mask);
        if isempty(t), continue; end

        [t,i] = sort(t);
        x   = dat.x(mask); x = x(i);
        via = dat.via(mask); via = via(i)/100;
        glc = dat.glc(mask); glc = glc(i);
        gln = dat.gln(mask); gln = gln(i);
        lac = dat.lac(mask); lac = lac(i);
        nh4 = dat.nh4(mask); nh4 = nh4(i);

        b = struct();
        b.name = sheets{s};
        b.t = t; b.x = x; b.via = via;
        b.glc = glc; b.gln = gln;
        b.lac = lac; b.nh4 = nh4;

        B{end+1} = b;
    end

    batches = B;

    set(hStatus,'String',...
        sprintf('%d개 배치 로드 완료',numel(batches)));
end

%% ============================================================
% 2) RUN MCMC (Adaptive)
% ============================================================
function run_mcmc(~,~)

    if isempty(batches)
        errordlg('엑셀 먼저 불러오세요.');
        return;
    end

    set(hRun,'Enable','off','BackgroundColor',[1 0 0],...
        'String','Running...'); drawnow;

    nIter = str2double(get(hIter,'String'));
    if isnan(nIter) || nIter < 500, nIter = 2000; end

    burn_frac = 0.5;
    burn = round(nIter*burn_frac);

    % -------- stoich fitting ----------
    S = fit_stoich(batches);

    % -------- MCMC 12개 초기값 ----------
    p = [ ...
        5, 5, ...            % Kglc, Kgln
        40, 40, ...          % KIamm, KIlac
        30, 30, ...          % KDamm, KDlac
        0.03, ...            % mglc
        1, 5, ...            % a1, a2
        0.01, ...            % dgln
        0.05, ...            % ramm
        0.1];                % QB1

    prop = 0.30 * p;       % bigger step to escape local region
    chain = zeros(nIter,numel(p));

    log_cur = total_loglik_full(p,S,batches);

    hW = waitbar(0,'Adaptive MCMC 진행중...');
    accept = 0;

    for i = 1:nIter
        p_prop = p + prop.*randn(size(p));

        if all(p_prop>0)
            log_p = total_loglik_full(p_prop,S,batches);
            if log(rand) < (log_p - log_cur)
                p = p_prop;
                log_cur = log_p;
                accept = accept + 1;
            end
        end

        chain(i,:) = p;

        % ---- Adaptive tuning ----
        if i < burn
            acc_rate = accept / i;
            if acc_rate > 0.40
                prop = prop * 1.05;
            elseif acc_rate < 0.20
                prop = prop * 0.95;
            end
        end

        if mod(i,20)==0
            waitbar(i/nIter,hW,...
                sprintf('Iter %d / %d (acc=%.2f)',i,nIter,accept/i));
        end
    end
    close(hW);

    posterior = chain(burn+1:end,:);
    p50 = median(posterior,1);

    plot_results(batches,S,p50);

    set(hRun,'Enable','on','BackgroundColor',[0.2 0.6 1],...
        'String','Run MCMC');
end

%% ============================================================
% Stoichiometric 6개 자동 추정
% ============================================================
function S = fit_stoich(B)

    tmax = max(cellfun(@(b)max(b.t),B));
    tG = linspace(0,tmax,200);

    fields = {'x','glc','gln','lac','nh4'};
    M = struct();

    for fi = 1:numel(fields)
        allInterp = nan(numel(B),numel(tG));
        for bb = 1:numel(B)
            t = B{bb}.t;
            v = B{bb}.(fields{fi});
            mask = ~isnan(t) & ~isnan(v);
            if sum(mask)>=2
                allInterp(bb,:) = interp1(t(mask),v(mask),tG,'linear','extrap');
            end
        end
        M.(fields{fi}) = mean(allInterp,1,'omitnan');
    end

    X = M.x; GLC = M.glc; GLN = M.gln; LAC = M.lac; NH4 = M.nh4;

    % umax
    idx = find(X > 0.2*max(X));
    idx = idx(1:round(0.4*numel(idx)));
    pv = polyfit(tG(idx),log(X(idx)),1);
    umax = pv(1);

    % kd
    idx2 = find(tG > 0.7*tmax);
    pd = polyfit(tG(idx2),log(X(idx2)),1);
    kd = -pd(1);

    % yields
    dX = max(X)-min(X);
    dGlc = min(GLC)-max(GLC);
    dGln = min(GLN)-max(GLN);
    dLac = max(LAC)-min(LAC);
    dNH4 = max(NH4)-min(NH4);

    Yxv_glc = dX/(abs(dGlc)+eps);
    Yxv_gln = dX/(abs(dGln)+eps);
    Ylac_glc = dLac/(abs(dGlc)+eps);
    Yamm_gln = dNH4/(abs(dGln)+eps);

    S = struct('umax',umax,'kd',kd,...
        'Yxv_glc',Yxv_glc,'Yxv_gln',Yxv_gln,...
        'Ylac_glc',Ylac_glc,'Yamm_gln',Yamm_gln);
end

%% ============================================================
% Total Log-Likelihood (weighted)
% ============================================================
function L = total_loglik_full(p,S,B)

    L = 0;

    for bb = 1:numel(B)
        b = B{bb};

        y0 = [b.x(1), b.glc(1), b.gln(1), b.lac(1), b.nh4(1), 1];

        try
            [t,y] = ode45(@(t,y) full_ode(t,y,p,S), b.t, y0);
        catch
            L = L - 1e8; return;
        end

        y = real(y);

        sigma_x  = 0.05 * range(b.x);
        sigma_glc = 0.05 * range(b.glc);
        sigma_gln = 0.05 * range(b.gln);
        sigma_lac = 0.05 * range(b.lac);
        sigma_nh4 = 0.05 * range(b.nh4);

        L = L - err(b.x, y(:,1), sigma_x);
        L = L - err(b.glc,y(:,2), sigma_glc);
        L = L - err(b.gln,y(:,3), sigma_gln);
        L = L - err(b.lac,y(:,4), sigma_lac);
        L = L - err(b.nh4,y(:,5), sigma_nh4);
    end
end

function LL = err(obs,sim,sigma)
    m = ~isnan(obs);
    if sum(m)<2, LL = 0; return; end
    r = (obs(m)-sim(m))/sigma;
    LL = 0.5*sum(r.^2) + sum(m)*log(sigma*sqrt(2*pi));
end

%% ============================================================
% FULL 6-ODE (기반식)
% ============================================================
function dy = full_ode(t,y,p,S)

    Xv = y(1); G=y(2); N=y(3);
    L=y(4); A=y(5); B1=y(6);

    % MCMC params
    Kglc=p(1); Kgln=p(2);
    KIamm=p(3); KIlac=p(4);
    KDamm=p(5); KDlac=p(6);
    mglc=p(7);
    a1=p(8); a2=p(9);
    dgln=p(10);
    ramm=p(11);
    QB1=p(12);

    % Stoich params
    umax=S.umax; kd=S.kd;
    Yxg=S.Yxv_glc;
    Yxn=S.Yxv_gln;
    Ylg=S.Ylac_glc;
    Yam=S.Yamm_gln;

    mu = umax * ...
        (G/(Kglc+G)) * ...
        (N/(Kgln+N)) * ...
        (KIlac/(KIlac+L)) * ...
        (KIamm/(KIamm+A));

    mu_d = kd * ...
        (L/(KDlac+L)) * ...
        (A/(KDamm+A));

    mgln = a1*N/(a2+N);

    dy = zeros(6,1);

    dy(1) = (mu-mu_d)*Xv;
    dy(2) = -((mu-mu_d)/Yxg + mglc)*Xv;
    dy(3) = -((mu-mu_d)/Yxn + mgln)*Xv - dgln*N;
    dy(4) = Ylg*((mu-mu_d)/Yxg + mglc)*Xv;
    dy(5) = Yam*((mu-mu_d)/Yxn)*Xv - ramm*Xv + dgln*N;
    dy(6) = QB1*Xv*(1-mu/umax)*B1;
end

%% ============================================================
% 결과 Plot (Batch별 곡선 = 각 배치 선 표시)
% ============================================================
function plot_results(B,S,p50)

figure('Name','Simulation Results','Color','w',...
    'Position',[1100 80 900 700]);

fields = {'x','glc','gln','lac','nh4'};
names = {'VCD','Glucose','Glutamine','Lactate','Ammonia'};
yl = {'10^6 cells/mL','mM','mM','mM','mM'};

cmap = lines(numel(B));

for i = 1:5
    subplot(3,2,i); hold on; grid on;

    for bb = 1:numel(B)
        b = B{bb};
        y0 = [b.x(1), b.glc(1), b.gln(1), b.lac(1), b.nh4(1), 1];
        [t,y] = ode45(@(t,y) full_ode(t,y,p50,S), b.t, y0);

        % batch별 파란선 → 각 배치 고유 컬러
        plot(b.t,y(:,i),'Color',cmap(bb,:),'LineWidth',1.4);
        plot(b.t,b.(fields{i}),'ko','MarkerSize',4);
    end

    xlabel('Time (h)');
    ylabel(yl{i});
    title(names{i});
end

% 파라미터 출력
subplot(3,2,6);
axis off;

txt = sprintf([ ...
    'Stoich (6):\n',...
    'umax=%.4f\nkd=%.4f\nYxv/glc=%.4f\nYxv/gln=%.4f\nYlac/glc=%.4f\nYamm/gln=%.4f\n\n',...
    'MCMC (12):\n',...
    'Kglc=%.3f Kgln=%.3f\nKIamm=%.3f KIlac=%.3f\nKDamm=%.3f KDlac=%.3f\n',...
    'mglc=%.3f a1=%.3f a2=%.3f\ndgln=%.3f ramm=%.3f QB1=%.3f\n'],...
    S.umax,S.kd,S.Yxv_glc,S.Yxv_gln,S.Ylac_glc,S.Yamm_gln,...
    p50(1),p50(2),p50(3),p50(4),p50(5),p50(6),...
    p50(7),p50(8),p50(9),p50(10),p50(11),p50(12));

text(0,0.5,txt,'FontSize',11,'VerticalAlignment','top');
end

%% ============================================================
% Utilities
% ============================================================
function v = toNum(c)
    if isnumeric(c), v=c; return; end
    v = str2double(string(c));
end

function dat = unitConvert(d,u)

    if contains(lower(u.t),'min')
        d.t = d.t/60;
    elseif contains(lower(u.t),'day')
        d.t = d.t*24;
    end

    MW.glc = 180.156;
    MW.gln = 146.14;
    MW.lac = 89.07;

    if contains(lower(u.glc),'g/')
        d.glc = d.glc*1000/MW.glc;
    end
    if contains(lower(u.gln),'g/')
        d.gln = d.gln*1000/MW.gln;
    end
    if contains(lower(u.lac),'g/')
        d.lac = d.lac*1000/MW.lac;
    end
    if contains(lower(u.nh4),'mg')||contains(lower(u.nh4),'nh4-n')
        d.nh4 = d.nh4/14.0067;
    end

    dat = d;
end

end


