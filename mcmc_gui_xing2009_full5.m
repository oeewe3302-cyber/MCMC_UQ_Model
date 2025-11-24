function plot_results(B,S,p50)

figure('Name','Simulation + Residual + d/dt Debug','Color','w',...
    'Position',[50 50 1500 900]);

fields = {'x','glc','gln','lac','nh4'};
names = {'VCD','Glucose','Glutamine','Lactate','Ammonia'};
yl = {'10^6 cells/mL','mM','mM','mM','mM'};

cmap = lines(numel(B));

%% --------------------------------------------
% 1) Batch별 실제값 + 예측값
% --------------------------------------------
for i = 1:5
    subplot(4,3,i); hold on; grid on;

    for bb = 1:numel(B)
        b = B{bb};
        y0 = [b.x(1), b.glc(1), b.gln(1), b.lac(1), b.nh4(1), 1];

        [t,y] = ode45(@(t,y) full_ode(t,y,p50,S), b.t, y0);

        plot(b.t,y(:,i),'Color',cmap(bb,:),'LineWidth',1.6);
        plot(b.t,b.(fields{i}),'ko','MarkerSize',4,'LineWidth',1.0);
    end
    xlabel('Time (h)'); ylabel(yl{i}); title(names{i});
end

%% --------------------------------------------
% 2) Residual Plot (obs - sim)
% --------------------------------------------
for i = 1:5
    subplot(4,3,i+5); hold on; grid on;

    for bb = 1:numel(B)
        b = B{bb};
        y0 = [b.x(1), b.glc(1), b.gln(1), b.lac(1), b.nh4(1), 1];
        [t,y] = ode45(@(t,y) full_ode(t,y,p50,S), b.t, y0);

        res = b.(fields{i}) - y(:,i);
        plot(b.t,res,'-o','Color',cmap(bb,:));
    end
    xlabel('Time (h)'); ylabel('Residual'); title([names{i} ' Residual']);
end

%% --------------------------------------------
% 3) d/dt 관계식 플롯 (Xing ODE 내부식)
% --------------------------------------------
% 예: 첫 번째 배치 기준으로 미분값 계산
b = B{1};
y0 = [b.x(1), b.glc(1), b.gln(1), b.lac(1), b.nh4(1), 1];
[t,y] = ode45(@(t,y) full_ode(t,y,p50,S), b.t, y0);

d = zeros(length(t),5);
for k = 1:length(t)
    dy = full_ode(t(k),y(k,:).',p50,S);
    d(k,:) = dy(1:5).';
end

dNames = {'dX/dt','dG/dt','dN/dt','dL/dt','dA/dt'};

for i = 1:5
    subplot(4,3,i+10); hold on; grid on;
    plot(t,d(:,i),'r-','LineWidth',1.5);
    xlabel('Time (h)'); ylabel(dNames{i});
    title(['Relation: ' dNames{i}]);
end

%% --------------------------------------------
% 4) 파라미터 텍스트
% --------------------------------------------
subplot(4,3,12);
axis off;

txt = sprintf([ ...
    'Stoich (6):\n',...
    'umax = %.4f\nkd = %.4f\nYx/glc = %.4f\nYx/gln = %.4f\nYlac/glc = %.4f\nYamm/gln = %.4f\n\n',...
    'MCMC (12):\n',...
    'Kglc = %.3f  Kgln = %.3f\nKIamm = %.3f  KIlac = %.3f\nKDamm = %.3f  KDlac = %.3f\n',...
    'mglc = %.3f  a1 = %.3f  a2 = %.3f\ndgln = %.3f  ramm = %.3f  QB1 = %.3f\n'],...
    S.umax,S.kd,S.Yxv_glc,S.Yxv_gln,S.Ylac_glc,S.Yamm_gln,...
    p50(1),p50(2),p50(3),p50(4),p50(5),p50(6),...
    p50(7),p50(8),p50(9),p50(10),p50(11),p50(12));

text(0,1,txt,'FontSize',10,'VerticalAlignment','top');

end
