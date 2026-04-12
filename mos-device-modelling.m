%% ========================================================================
%  MOSFET SIMULATION: Square Law vs BSIM Model Under Short-Channel Effects
%  Author  : Srilakshmi Badri Seshadri (230907134)
%  Project : Breakdown of Square Law MOS Model vs BSIM Model


clear; clc; close all;


%  GLOBAL DEVICE PARAMETERS (BSIM-inspired, 65 nm process node)

q       = 1.602e-19;        % Elementary charge [C]
eps0    = 8.854e-12;        % Permittivity of free space [F/m]
eps_si  = 11.7 * eps0;      % Silicon permittivity [F/m]
eps_ox  = 3.9  * eps0;      % SiO2 permittivity [F/m]
k_B     = 1.381e-23;        % Boltzmann constant [J/K]
T       = 300;              % Temperature [K]
Vt_th   = k_B * T / q;     % Thermal voltage ~0.02585 V

% Oxide & threshold
tox     = 2e-9;             % Gate oxide thickness [m] (2 nm, 65-nm node)
Cox     = eps_ox / tox;     % Oxide capacitance per unit area [F/m^2]
Vth0    = 0.45;             % Long-channel threshold voltage [V]
W       = 1e-6;             % Width [m]

% Mobility
mu_n0   = 0.04;             % Low-field electron mobility [m^2/V·s] (400 cm^2/Vs)
theta   = 1.2;              % Mobility degradation coefficient [1/V]

% Velocity saturation
vsat    = 1e5;              % Saturation velocity [m/s]  (1e7 cm/s)
Ecrit   = vsat / mu_n0;    % Critical electric field [V/m]

% BSIM short-channel parameters
Xj      = 20e-9;            % Junction depth [m]
lambda0 = 0.05;             % Output conductance parameter (long-channel)

% Subthreshold parameters
n_id    = 1.3;              % Ideality factor (long channel)
n_id_sc = 1.8;              % Ideality (short channel, DIBL degrades SS)

% DIBL parameters
DIBL_coeff = 0.35;          % Empirical DIBL coefficient [V/V]
L_ref      = 100e-9;        % Reference length for DIBL scaling [m]

% Channel lengths to simulate
L_vals  = [1000 500 250 100 65 45] * 1e-9;   % [m]
L_nm    = L_vals * 1e9;                        % [nm] for labels
colors  = lines(length(L_vals));

% Voltage sweeps
Vgs_range   = 0 : 0.01 : 1.2;   % [V]
Vds_lin     = 0 : 0.01 : 1.2;   % [V]
Vds_low     = 0.05;              % Low Vds for DIBL extraction
Vds_high    = 1.0;               % High Vds for DIBL extraction



%  5.1 — LONG-CHANNEL SQUARE LAW: Id-Vgs and Id-Vds

fprintf('==> Plotting Section 5.1: Square Law MOSFET\n');

L_long  = 1000e-9;   % 1 µm long channel
mu_n    = mu_n0;     % no mobility degradation for ideal model

figure('Name','5.1 Square Law - Id-Vgs','NumberTitle','off',...
       'Position',[50 50 700 500]);

% Id-Vgs (Vds = 1.0 V, saturation)
Vds_fixed = 1.0;
Id_sq = zeros(size(Vgs_range));
for i = 1:length(Vgs_range)
    Vgs = Vgs_range(i);
    Vov = Vgs - Vth0;
    if Vov <= 0
        Id_sq(i) = 0;
    elseif Vds_fixed >= Vov   % saturation
        Id_sq(i) = 0.5 * mu_n * Cox * (W/L_long) * Vov^2;
    else                       % linear
        Id_sq(i) = mu_n * Cox * (W/L_long) * ((Vov)*Vds_fixed - 0.5*Vds_fixed^2);
    end
end

subplot(1,2,1);
plot(Vgs_range, Id_sq*1e6, 'b-', 'LineWidth', 2);
xlabel('V_{GS} (V)','FontSize',12);
ylabel('I_D (\muA)','FontSize',12);
title('Square Law: I_D vs V_{GS}','FontSize',13,'FontWeight','bold');
subtitle(sprintf('L = %d nm, V_{DS} = %.1f V', L_long*1e9, Vds_fixed));
grid on; box on;
xline(Vth0,'r--','V_{th}','LabelVerticalAlignment','bottom','FontSize',10);
legend('I_D (Square Law)','Location','northwest');

% Id-Vds family (square law)
subplot(1,2,2);
Vgs_family = [0.4 0.6 0.8 1.0 1.2];
hold on;
for v = 1:length(Vgs_family)
    Vgs   = Vgs_family(v);
    Vov   = Vgs - Vth0;
    Id_vds = zeros(size(Vds_lin));
    for i = 1:length(Vds_lin)
        Vds = Vds_lin(i);
        if Vov <= 0
            Id_vds(i) = 0;
        elseif Vds < Vov
            Id_vds(i) = mu_n*Cox*(W/L_long)*((Vov)*Vds - 0.5*Vds^2);
        else
            Id_vds(i) = 0.5*mu_n*Cox*(W/L_long)*Vov^2*(1+lambda0*(Vds-Vov));
        end
    end
    plot(Vds_lin, Id_vds*1e6, 'LineWidth', 2, 'DisplayName', ...
         sprintf('V_{GS} = %.1f V', Vgs));
end
xlabel('V_{DS} (V)','FontSize',12);
ylabel('I_D (\muA)','FontSize',12);
title('Square Law: I_D vs V_{DS}','FontSize',13,'FontWeight','bold');
subtitle(sprintf('L = %d nm', L_long*1e9));
legend('Location','southeast','FontSize',9);
grid on; box on;

sgtitle('Section 5.1: Long-Channel Square-Law MOSFET','FontSize',14,...
        'FontWeight','bold');

%  5.3 — Id-Vgs CHARACTERISTICS (Multiple channel lengths, BSIM-style)

fprintf('==> Plotting Section 5.3: Id-Vgs Multi-Length\n');

figure('Name','5.3 Id-Vgs Multi-Length','NumberTitle','off',...
       'Position',[50 600 800 500]);
hold on;

Id_linear_all = zeros(length(L_vals), length(Vgs_range));  % store for later

for k = 1:length(L_vals)
    L      = L_vals(k);
    Vth_k  = Vth_bsim(Vth0, L, Xj, Vds_low, DIBL_coeff, L_ref);
    Id_arr = zeros(size(Vgs_range));
    for i = 1:length(Vgs_range)
        Id_arr(i) = Id_bsim(Vgs_range(i), Vds_high, Vth_k, L, W, ...
                            mu_n0, Cox, theta, vsat, lambda0);
    end
    Id_linear_all(k,:) = Id_arr;
    plot(Vgs_range, Id_arr*1e6, 'Color', colors(k,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('L = %d nm', round(L_nm(k))));
end

% Overlay square law for comparison
plot(Vgs_range, Id_sq*1e6, 'k--', 'LineWidth', 1.5, ...
     'DisplayName', 'Square Law (L=1µm)');

xlabel('V_{GS} (V)','FontSize',12);
ylabel('I_D (\muA)','FontSize',12);
title('Section 5.3: I_D–V_{GS} Characteristics vs Channel Length','FontSize',13,...
      'FontWeight','bold');
subtitle(sprintf('V_{DS} = %.1f V, W = 1 µm', Vds_high));
legend('Location','northwest','FontSize',9);
grid on; box on;
xlim([0 1.2]); ylim([0 inf]);


%  5.4 — LOG(Id)-Vgs PLOT & SUBTHRESHOLD SLOPE EXTRACTION

fprintf('==> Plotting Section 5.4: Log(Id)-Vgs & Subthreshold Slope\n');

figure('Name','5.4 Log(Id)-Vgs','NumberTitle','off',...
       'Position',[850 50 800 500]);

% Build subthreshold + above-threshold current (BSIM-like with subVT)
Vgs_fine = -0.1 : 0.005 : 1.2;
hold on;

SS_vals = zeros(1, length(L_vals));   % Subthreshold slope per length

for k = 1:length(L_vals)
    L     = L_vals(k);
    % DIBL shifts Vth at Vds_high
    Vth_k = Vth_bsim(Vth0, L, Xj, Vds_high, DIBL_coeff, L_ref);
    % Effective ideality factor degrades with short channel
    n_eff = n_id + (n_id_sc - n_id) * (1 - min(L/L_vals(1),1));
    SS_vals(k) = n_eff * Vt_th * log(10) * 1000;  % mV/dec

    Id_log = zeros(size(Vgs_fine));
    for i = 1:length(Vgs_fine)
        Vgs = Vgs_fine(i);
        % Subthreshold component
        Id_sub = (W/L) * mu_n0 * Cox * (n_eff-1) * Vt_th^2 * ...
                  exp((Vgs - Vth_k)/(n_eff * Vt_th)) * ...
                  (1 - exp(-Vds_high / Vt_th));
        % Above-threshold component
        Id_above = Id_bsim(Vgs, Vds_high, Vth_k, L, W, ...
                           mu_n0, Cox, theta, vsat, lambda0);
        Id_log(i) = Id_sub + Id_above;
    end
    Id_log = max(Id_log, 1e-18);  % floor for log
    semilogy(Vgs_fine, Id_log, 'Color', colors(k,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('L=%d nm, SS≈%.0f mV/dec', ...
             round(L_nm(k)), SS_vals(k)));
end

xlabel('V_{GS} (V)','FontSize',12);
ylabel('I_D (A)','FontSize',12);
title('Section 5.4: log(I_D)–V_{GS}: Subthreshold Region','FontSize',13,...
      'FontWeight','bold');
subtitle(sprintf('V_{DS} = %.1f V — SS degrades with L scaling', Vds_high));
legend('Location','southeast','FontSize',9);
grid on; box on;
ylim([1e-16 1e-3]);

% ---- Subthreshold Slope bar chart ----
figure('Name','5.4 Subthreshold Slope vs L','NumberTitle','off',...
       'Position',[850 600 600 400]);
bar(1:length(L_vals), SS_vals, 'FaceColor',[0.2 0.6 0.9],'EdgeColor','k');
yline(60,'r--','Ideal Limit (60 mV/dec)','LabelHorizontalAlignment','left',...
      'FontSize',10,'LineWidth',1.5);
set(gca,'XTickLabel', arrayfun(@(x) sprintf('%d nm',x), round(L_nm),...
    'UniformOutput',false));
xlabel('Channel Length','FontSize',12);
ylabel('Subthreshold Slope (mV/dec)','FontSize',12);
title('Section 5.4: Subthreshold Slope vs Channel Length','FontSize',13,...
      'FontWeight','bold');
grid on; ylim([50 120]);


%  5.5 — VELOCITY SATURATION (Id-Vds & Id-Vgs: short vs long channel)

fprintf('==> Plotting Section 5.5: Velocity Saturation\n');

figure('Name','5.5 Velocity Saturation','NumberTitle','off',...
       'Position',[50 50 1100 500]);


subplot(1,2,1);
L_long_vs  = 1000e-9;
L_short_vs = 65e-9;
Vth_long   = Vth_bsim(Vth0, L_long_vs,  Xj, Vds_high, DIBL_coeff, L_ref);
Vth_short  = Vth_bsim(Vth0, L_short_vs, Xj, Vds_high, DIBL_coeff, L_ref);

Id_long_vgs  = zeros(size(Vgs_range));
Id_short_vgs = zeros(size(Vgs_range));
Id_sq_vsat   = zeros(size(Vgs_range));   % square law (no vsat)

for i = 1:length(Vgs_range)
    Vgs = Vgs_range(i);
    Id_long_vgs(i)  = Id_bsim(Vgs, Vds_high, Vth_long,  L_long_vs,  W,...
                               mu_n0, Cox, theta, vsat, lambda0);
    Id_short_vgs(i) = Id_bsim(Vgs, Vds_high, Vth_short, L_short_vs, W,...
                               mu_n0, Cox, theta, vsat, lambda0);
    % Square law (no velocity sat)
    Vov = Vgs - Vth0;
    if Vov > 0
        Id_sq_vsat(i) = 0.5*mu_n0*Cox*(W/L_long_vs)*Vov^2;
    end
end

hold on;
plot(Vgs_range, Id_sq_vsat*1e6,  'k--', 'LineWidth',1.5,'DisplayName','Square Law (L=1µm)');
plot(Vgs_range, Id_long_vgs*1e6, 'b-',  'LineWidth',2,  'DisplayName','BSIM L=1000 nm');
plot(Vgs_range, Id_short_vgs*1e6,'r-',  'LineWidth',2,  'DisplayName','BSIM L=65 nm');

xlabel('V_{GS} (V)','FontSize',12);
ylabel('I_D (\muA)','FontSize',12);
title('Velocity Saturation: I_D vs V_{GS}','FontSize',12,'FontWeight','bold');
subtitle('Quadratic→Linear transition in short-channel');
legend('Location','northwest','FontSize',9); grid on; box on;

% Annotate exponent change
Vov_annot = 0.9 - Vth_short;
if Vov_annot > 0
    text(0.85, Id_short_vgs(end)*1e6*0.7, '\alpha \approx 1 (velocity sat.)',...
         'Color','r','FontSize',10,'FontAngle','italic');
    text(0.85, Id_long_vgs(end)*1e6*0.6,  '\alpha \approx 2 (square law)',...
         'Color','b','FontSize',10,'FontAngle','italic');
end


subplot(1,2,2);
hold on;
Vgs_test = 1.0;
for k = [1 3 5 6]   % Select a few lengths
    L     = L_vals(k);
    Vth_k = Vth_bsim(Vth0, L, Xj, 0, DIBL_coeff, L_ref);
    Id_vd = zeros(size(Vds_lin));
    for i = 1:length(Vds_lin)
        Id_vd(i) = Id_bsim(Vgs_test, Vds_lin(i), Vth_k, L, W, ...
                            mu_n0, Cox, theta, vsat, lambda0);
    end
    plot(Vds_lin, Id_vd*1e6, 'Color', colors(k,:), 'LineWidth', 2,...
         'DisplayName', sprintf('L = %d nm', round(L_nm(k))));
end
% Square law Id-Vds for reference
Vov_ref = Vgs_test - Vth0;
Id_sq_vds = zeros(size(Vds_lin));
for i = 1:length(Vds_lin)
    Vds = Vds_lin(i);
    if Vov_ref <= 0
        Id_sq_vds(i) = 0;
    elseif Vds < Vov_ref
        Id_sq_vds(i) = mu_n0*Cox*(W/L_long_vs)*((Vov_ref)*Vds - 0.5*Vds^2);
    else
        Id_sq_vds(i) = 0.5*mu_n0*Cox*(W/L_long_vs)*Vov_ref^2*(1+lambda0*(Vds-Vov_ref));
    end
end
plot(Vds_lin, Id_sq_vds*1e6, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Square Law (L=1µm)');

xlabel('V_{DS} (V)','FontSize',12);
ylabel('I_D (\muA)','FontSize',12);
title('Velocity Saturation: I_D vs V_{DS}','FontSize',12,'FontWeight','bold');
subtitle(sprintf('V_{GS} = %.1f V — saturation occurs at lower V_{DS} for short L', Vgs_test));
legend('Location','southeast','FontSize',9); grid on; box on;

sgtitle('Section 5.5: Velocity Saturation Effects','FontSize',14,'FontWeight','bold');


%  5.2 & 5.6 — Vth vs CHANNEL LENGTH (Threshold Voltage Roll-off)

fprintf('==> Plotting Section 5.2 & 5.6: Vth Roll-off vs L\n');

L_sweep  = logspace(-8, -6, 200);   % 10 nm to 1000 nm
Vth_low  = arrayfun(@(L) Vth_bsim(Vth0, L, Xj, Vds_low,  DIBL_coeff, L_ref), L_sweep);
Vth_high = arrayfun(@(L) Vth_bsim(Vth0, L, Xj, Vds_high, DIBL_coeff, L_ref), L_sweep);
Vth_sq   = Vth0 * ones(size(L_sweep));  % Square law: constant

figure('Name','5.6 Vth Roll-off','NumberTitle','off',...
       'Position',[50 50 750 520]);
semilogx(L_sweep*1e9, Vth_sq,   'k--',  'LineWidth', 1.5, 'DisplayName','Square Law (constant)');
hold on;
semilogx(L_sweep*1e9, Vth_low,  'b-',   'LineWidth', 2,   'DisplayName', ...
         sprintf('BSIM: V_{DS}=%.2f V (low)', Vds_low));
semilogx(L_sweep*1e9, Vth_high, 'r-',   'LineWidth', 2,   'DisplayName', ...
         sprintf('BSIM: V_{DS}=%.1f V (high, DIBL)', Vds_high));

% Mark the simulated L points
for k = 1:length(L_vals)
    Vth_pt_lo = Vth_bsim(Vth0, L_vals(k), Xj, Vds_low,  DIBL_coeff, L_ref);
    Vth_pt_hi = Vth_bsim(Vth0, L_vals(k), Xj, Vds_high, DIBL_coeff, L_ref);
    semilogx(L_nm(k), Vth_pt_lo, 'bs', 'MarkerSize', 8, 'MarkerFaceColor','b','HandleVisibility','off');
    semilogx(L_nm(k), Vth_pt_hi, 'r^', 'MarkerSize', 8, 'MarkerFaceColor','r','HandleVisibility','off');
end

xlabel('Channel Length L (nm)','FontSize',12);
ylabel('Threshold Voltage V_{th} (V)','FontSize',12);
title('Sections 5.2 & 5.6: Threshold Voltage Roll-off vs Channel Length','FontSize',13,...
      'FontWeight','bold');
subtitle('Symbols = simulated nodes; curves = analytical BSIM model');
legend('Location','southwest','FontSize',10);
grid on; box on;
xlim([10 1100]);
ylim([0 0.6]);

% Annotate roll-off region



%  5.7 — DIBL ANALYSIS (DIBL vs Channel Length)

fprintf('==> Plotting Section 5.7: DIBL Analysis\n');

Vth_at_low  = arrayfun(@(L) Vth_bsim(Vth0, L, Xj, Vds_low,  DIBL_coeff, L_ref), L_sweep);
Vth_at_high = arrayfun(@(L) Vth_bsim(Vth0, L, Xj, Vds_high, DIBL_coeff, L_ref), L_sweep);
DIBL_sweep  = (Vth_at_low - Vth_at_high) / (Vds_high - Vds_low) * 1000;  % mV/V

DIBL_nodes  = zeros(1, length(L_vals));
for k = 1:length(L_vals)
    Vt_lo = Vth_bsim(Vth0, L_vals(k), Xj, Vds_low,  DIBL_coeff, L_ref);
    Vt_hi = Vth_bsim(Vth0, L_vals(k), Xj, Vds_high, DIBL_coeff, L_ref);
    DIBL_nodes(k) = (Vt_lo - Vt_hi) / (Vds_high - Vds_low) * 1000;
end

figure('Name','5.7 DIBL Analysis','NumberTitle','off',...
       'Position',[850 50 1100 500]);


subplot(1,2,1);
semilogx(L_sweep*1e9, DIBL_sweep, 'b-', 'LineWidth', 2, 'DisplayName','DIBL (mV/V)');
hold on;
semilogx(L_nm, DIBL_nodes, 'ro', 'MarkerSize', 9, 'MarkerFaceColor','r',...
         'DisplayName','Simulated nodes');
yline(100,'g--','100 mV/V limit','LabelHorizontalAlignment','left',...
      'FontSize',10,'LineWidth',1.5);
yline(200,'r--','200 mV/V (45 nm limit)','LabelHorizontalAlignment','left',...
      'FontSize',10,'LineWidth',1.5);

xlabel('Channel Length L (nm)','FontSize',12);
ylabel('DIBL (mV/V)','FontSize',12);
title('Section 5.7: DIBL vs Channel Length','FontSize',13,'FontWeight','bold');
subtitle(sprintf('DIBL = [V_{th}(V_{DS,lo}) - V_{th}(V_{DS,hi})] / \\DeltaV_{DS}'));
legend('Location','northeast','FontSize',9);
grid on; box on;


subplot(1,2,2);
Vth_lo_nodes = arrayfun(@(L) Vth_bsim(Vth0,L,Xj,Vds_low, DIBL_coeff,L_ref), L_vals);
Vth_hi_nodes = arrayfun(@(L) Vth_bsim(Vth0,L,Xj,Vds_high,DIBL_coeff,L_ref), L_vals);

x_bar = 1:length(L_vals);
bar(x_bar, [Vth_lo_nodes; Vth_hi_nodes]'*1000, 'grouped');
set(gca,'XTickLabel', arrayfun(@(x) sprintf('%d nm',x), round(L_nm),...
    'UniformOutput',false));
legend(sprintf('V_{DS} = %.2f V', Vds_low), sprintf('V_{DS} = %.1f V', Vds_high),...
       'Location','northeast','FontSize',9);
xlabel('Channel Length','FontSize',12);
ylabel('V_{th} (mV)','FontSize',12);
title('Section 5.7: V_{th} Shift Due to DIBL','FontSize',13,'FontWeight','bold');
subtitle('DIBL lowers V_{th} more severely at short L and high V_{DS}');
grid on; box on;


fprintf('==> Plotting Summary: Square Law vs BSIM\n');

figure('Name','Summary - Square Law vs BSIM','NumberTitle','off',...
       'Position',[50 50 1300 700]);

L_compare = 65e-9;
Vth_bsim65 = Vth_bsim(Vth0, L_compare, Xj, Vds_high, DIBL_coeff, L_ref);


Id_sq_comp = zeros(size(Vgs_range));
Id_bsim_comp = zeros(size(Vgs_range));
for i = 1:length(Vgs_range)
    Vgs = Vgs_range(i);
    Vov = Vgs - Vth0;
    if Vov > 0
        Id_sq_comp(i) = 0.5*mu_n0*Cox*(W/L_long)*Vov^2;
    end
    Id_bsim_comp(i) = Id_bsim(Vgs, Vds_high, Vth_bsim65, L_compare, W,...
                               mu_n0, Cox, theta, vsat, lambda0);
end

% --- Subplot 1: Id-Vgs ---
subplot(2,3,1);
plot(Vgs_range, Id_sq_comp*1e6,  'b--', 'LineWidth',2,'DisplayName','Square Law');
hold on;
plot(Vgs_range, Id_bsim_comp*1e6,'r-',  'LineWidth',2,'DisplayName','BSIM 65nm');
xlabel('V_{GS} (V)'); ylabel('I_D (\muA)');
title('I_D–V_{GS}: Linear Scale'); legend; grid on;
xline(Vth0,'b:','V_{th,SL}','FontSize',8);
xline(Vth_bsim65,'r:','V_{th,BSIM}','FontSize',8);

% --- Subplot 2: Log Id-Vgs ---
subplot(2,3,2);
n_eff65 = n_id + (n_id_sc-n_id)*(1-min(L_compare/L_vals(1),1));
Id_log65 = zeros(size(Vgs_fine));
for i = 1:length(Vgs_fine)
    Vgs = Vgs_fine(i);
    Id_sub = (W/L_compare)*mu_n0*Cox*(n_eff65-1)*Vt_th^2 * ...
              exp((Vgs-Vth_bsim65)/(n_eff65*Vt_th))*(1-exp(-Vds_high/Vt_th));
    Id_above = Id_bsim(Vgs, Vds_high, Vth_bsim65, L_compare, W,...
                       mu_n0, Cox, theta, vsat, lambda0);
    Id_log65(i) = max(Id_sub + Id_above, 1e-18);
end
Id_sq_log = zeros(size(Vgs_fine));
for i = 1:length(Vgs_fine)
    Vov = Vgs_fine(i)-Vth0;
    if Vov>0; Id_sq_log(i)=0.5*mu_n0*Cox*(W/L_long)*Vov^2; end
    Id_sq_log(i) = max(Id_sq_log(i),1e-18);
end
semilogy(Vgs_fine, Id_sq_log,'b--','LineWidth',2,'DisplayName','Square Law');
hold on;
semilogy(Vgs_fine, Id_log65,'r-','LineWidth',2,'DisplayName','BSIM 65nm');
xlabel('V_{GS} (V)'); ylabel('I_D (A)');
title('log(I_D)–V_{GS}: Subthreshold'); legend; grid on;
ylim([1e-16 1e-3]);

% --- Subplot 3: Id-Vds ---
subplot(2,3,3);
hold on;
for v = [2 4 5]
    Vgs_t = Vgs_family(v);
    Vov = Vgs_t - Vth0;
    Vth65_vds = Vth_bsim(Vth0, L_compare, Xj, 0, DIBL_coeff, L_ref);
    Id_s  = zeros(size(Vds_lin));
    Id_b  = zeros(size(Vds_lin));
    for i = 1:length(Vds_lin)
        if Vov<=0; Id_s(i)=0;
        elseif Vds_lin(i)<Vov
            Id_s(i)=mu_n0*Cox*(W/L_long)*((Vov)*Vds_lin(i)-0.5*Vds_lin(i)^2);
        else
            Id_s(i)=0.5*mu_n0*Cox*(W/L_long)*Vov^2*(1+lambda0*(Vds_lin(i)-Vov));
        end
        Id_b(i) = Id_bsim(Vgs_t, Vds_lin(i), Vth65_vds, L_compare, W,...
                          mu_n0, Cox, theta, vsat, lambda0);
    end
    plot(Vds_lin, Id_s*1e6,'--','LineWidth',1.5,...
         'DisplayName',sprintf('SL V_{GS}=%.1f',Vgs_t));
    plot(Vds_lin, Id_b*1e6,'-','LineWidth',1.5,...
         'DisplayName',sprintf('BSIM V_{GS}=%.1f',Vgs_t));
end
xlabel('V_{DS} (V)'); ylabel('I_D (\muA)');
title('I_D–V_{DS}: Square Law (--) vs BSIM (-)'); legend('FontSize',7); grid on;

% --- Subplot 4: Vth roll-off ---
subplot(2,3,4);
semilogx(L_sweep*1e9, Vth_sq,   'k--','LineWidth',1.5,'DisplayName','Square Law');
hold on;
semilogx(L_sweep*1e9, Vth_low,  'b-', 'LineWidth',2,'DisplayName','BSIM V_{DS,lo}');
semilogx(L_sweep*1e9, Vth_high, 'r-', 'LineWidth',2,'DisplayName','BSIM V_{DS,hi}');
xlabel('L (nm)'); ylabel('V_{th} (V)');
title('V_{th} Roll-off'); legend; grid on; xlim([10 1100]);

% --- Subplot 5: DIBL vs L ---
subplot(2,3,5);
semilogx(L_sweep*1e9, DIBL_sweep,'b-','LineWidth',2);
hold on;
semilogx(L_nm, DIBL_nodes,'ro','MarkerSize',8,'MarkerFaceColor','r');
yline(100,'g--','100 mV/V'); yline(200,'r--','200 mV/V');
xlabel('L (nm)'); ylabel('DIBL (mV/V)');
title('DIBL vs Channel Length'); grid on;

% --- Subplot 6: SS vs L ---
subplot(2,3,6);
bar(1:length(L_vals), SS_vals, 'FaceColor',[0.2 0.7 0.3]);
yline(60,'r--','Ideal 60 mV/dec','FontSize',9,'LineWidth',1.5);
set(gca,'XTickLabel',arrayfun(@(x) sprintf('%d',x),round(L_nm),'UniformOutput',false));
xlabel('Channel Length (nm)'); ylabel('SS (mV/dec)');
title('Subthreshold Slope vs L'); grid on;

sgtitle('SUMMARY: Square-Law Model vs BSIM — Short-Channel Effects',...
        'FontSize',15,'FontWeight','bold');


fprintf('\n=====================================================\n');
fprintf(' EXTRACTED DEVICE PARAMETERS TABLE\n');
fprintf('=====================================================\n');
fprintf('%-12s %-10s %-12s %-12s %-12s\n',...
        'L (nm)','Vth (V)','DIBL (mV/V)','SS (mV/dec)','Ion (µA)');
fprintf('-----------------------------------------------------\n');
for k = 1:length(L_vals)
    Vt_lo = Vth_bsim(Vth0, L_vals(k), Xj, Vds_low,  DIBL_coeff, L_ref);
    Vt_hi = Vth_bsim(Vth0, L_vals(k), Xj, Vds_high, DIBL_coeff, L_ref);
    D_val = (Vt_lo - Vt_hi)/(Vds_high-Vds_low)*1000;
    Ion   = Id_bsim(1.2, Vds_high, Vt_hi, L_vals(k), W, mu_n0, Cox,...
                    theta, vsat, lambda0)*1e6;
    fprintf('%-12.0f %-10.4f %-12.1f %-12.1f %-12.2f\n',...
            L_nm(k), Vt_hi, D_val, SS_vals(k), Ion);
end
fprintf('=====================================================\n\n');

fprintf('All simulations complete. %d figures generated.\n', length(findobj('Type','figure')));



function Vth = Vth_bsim(Vth0, L, Xj, Vds, DIBL_coeff, L_ref)

    delta_Vth_rolloff = 0.15 * (Xj / L) * (1 + Vds * 0.3);
    delta_Vth_rolloff = min(delta_Vth_rolloff, Vth0 * 0.85);  % cap

    % DIBL shift
    DIBL_factor = DIBL_coeff * (L_ref / L)^1.5;   % increases at short L
    DIBL_shift  = DIBL_factor * Vds;

    Vth = Vth0 - delta_Vth_rolloff - DIBL_shift;
    Vth = max(Vth, 0.01);  % physical floor
end

function Id = Id_bsim(Vgs, Vds, Vth, L, W, mu_n0, Cox, theta, vsat, lambda0)

    Vov = Vgs - Vth;
    if Vov <= 0
        Id = 0;
        return;
    end

    % Effective mobility with degradation
    mu_eff = mu_n0 / (1 + theta * Vov);

    % Velocity-saturation saturating Vds
    Ecrit  = vsat / mu_eff;
    Vdsat  = min(Vov, Ecrit * L);       % Vdsat limited by velocity sat
    Vds_eff = min(Vds, Vdsat);          % clamp Vds to Vdsat

    % Drain current (unified linear/saturation with velocity sat)
    Id_lin = mu_eff * Cox * (W/L) * (Vov * Vds_eff - 0.5 * Vds_eff^2) ...
             / (1 + Vds_eff / (Ecrit * L));

    % Output conductance (channel length modulation in saturation)
    if Vds > Vdsat
        lambda = lambda0 / (L * 1e6 + 0.01);  % lambda increases with shorter L
        lambda = min(lambda, 0.5);
        Id = Id_lin * (1 + lambda * (Vds - Vdsat));
    else
        Id = Id_lin;
    end

    Id = max(Id, 0);
end