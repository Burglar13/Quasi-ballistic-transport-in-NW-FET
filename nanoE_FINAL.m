%%% Quasi-Ballistic Transport in a nW-FET   %%%
%%% Tommaso Ugolini                         %%%
%%% University of Bologna                   %%%
%%% February 2023                           %%%

clear 
clc
close all


R = 3e-9;                       L_g = 20e-9;    % = z2 - z1
t_ox = 1e-9;                    N_D = 1e20; % 10^20 is in cm^-3
epsilon_s = 11.68 * 8.854e-12;  epsilon_ox = 3.9 * 8.854e-12;
k_B = physconst('Boltzmann');   T = 300;    % °K
q = 1.602176634e-19;            % electron charge [Coulomb]
E_g = q .* 1.12;                % 1.12 [eV] energy gap of silicon
n_i = 1.5e10;    % 9.65e9;      % [cm^-3] intrinsic carrier concentration
h = 6.62607015e-34;             h_bar = h/(2*pi);
m_0 = 9.1093837e-31;
m_t = 0.19 * m_0;               m_l = 0.98 * m_0;

z_1 = 0;    z_2 = L_g;
% E_0_eV = 0:0.025:0.2;      E_0 = q .* E_0_eV;  % from 0 to 0.2 [eV] 
lambda_0 = 38e-9;

z = (z_1 + 1e-9):1e-9:z_2;

% E_F = E_g/2 + k_B * T * log(N_D / n_i);
% E_F = k_B * T * log(N_D / n_i);
E_F = E_g/2;

% u_F = E_g / (2 * k_B * T);
u_F = E_F / (k_B * T);

v = 0:10:40;                    V_DS = v .* k_B .* T ./ q;
u_g = 0:3.75:15;                V_g = u_g .* k_B .* T ./ q;

E_FS = 0;       E_FD = k_B .* T .* v;               

epsilon_0 = 0:0.2:20; E_0 = epsilon_0 .* (k_B * T) ./ q;

m_r = (m_t + 2 * m_l * m_t) ./ 3;

gamma = (2 .* epsilon_s ./ epsilon_ox) .* log( 1 + (t_ox ./ R) );

lambda_c = R/2 * sqrt(1 + gamma);

lambda_D = sqrt( (epsilon_s * k_B * T) / (q.^2 .* (N_D .* 1e6)) );

lc_ld_2 = (lambda_c/lambda_D).^2;
ld_lc_2 = (lambda_D/lambda_c).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  uc1 - ug , uc2 - ug

szu = length(u_g);  szv = length(v);

uc1_ug = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        uc1_ug(ii,jj) = lc_ld_2 .* ( sqrt( 1 + 2 .* ld_lc_2 .* (u_F - u_g(ii)) ) - 1 );      
    end
end

uc1 = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        uc1(ii,jj) = uc1_ug(ii,jj) + u_g(ii);
    end
end
    
figure;
plot(u_g,uc1,'LineWidth',5);
title("$u_{c1}$ al variare di $u_g$");
xlabel("$u_g$");ylabel("$u_{c1}$");
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on


uc2_ug = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        uc2_ug(ii,jj) = lc_ld_2 .* ( sqrt( 1 + 2 .* ld_lc_2 .* (u_F + v(jj) - u_g(ii)) ) - 1 );
    end
end

uc2 = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        uc2(ii,jj) = uc2_ug(ii,jj) + u_g(ii);
    end
end

figure;
subplot(1,2,1);
plot(u_g,uc2,'LineWidth',5);
title("$u_{c2}$ al variare di $u_g$");
xlabel("$u_g$");ylabel("$u_{c2}$");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
subplot(1,2,2);
plot(v,uc2.','LineWidth',5);
title("$u_{c2}$ al variare di $v$");
xlabel("$v$");ylabel("$u_{c2}$");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on

S = uc2_ug ./ uc1_ug;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  z_m

z_m = z_1 + (lambda_c / 2) .* log( (exp(L_g ./ lambda_c) - S) ./ (S - exp(-L_g ./ lambda_c)) );

figure;
subplot(1,2,1);
plot(u_g,z_m,'LineWidth',5);
title("$z_{m}$ al variare di $u_g$");
xlabel("$u_g$");ylabel("$z_{m}$ [m]");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
subplot(1,2,2);
plot(v,z_m.','LineWidth',5);
title("$z_{m}$ al variare di $v$");
xlabel("$v$");ylabel("$z_{m}$ [m]");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% u_cm

u_cm = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        u_cm(ii,jj) = u_g(ii)  +  uc2_ug(ii,jj) .* ( sinh(z_m(ii,jj) ./ lambda_c) ./ sinh(L_g ./ lambda_c) )  + ...
            uc1_ug(ii,jj) .* ( sinh((L_g - z_m(ii,jj)) ./ lambda_c) ./ sinh(L_g ./ lambda_c) );
    end
end

figure;
subplot(1,2,1);
plot(u_g,u_cm,'LineWidth',5);
title("$u_{c}(z_m)$ al variare di $u_g$");
xlabel("$u_g$");ylabel("$u_{c}(z_m)$");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
subplot(1,2,2);
plot(v,u_cm.','LineWidth',5);
title("$u_{c}(z_m)$ al variare di $v$");
xlabel("$v$");ylabel("$u_{c}(z_m)$");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L_kT

L_kT = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        L_kT(ii,jj) = ( sqrt(2) .* lambda_c )  ./  sqrt( u_cm(ii,jj) - u_g(ii) );
    end
end

figure;
subplot(1,2,1);
plot(u_g,L_kT,'LineWidth',5);
title("$L_{kT}$ al variare di $u_g$");
xlabel("$u_g$");ylabel("$L_{kT}$ [m]");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
subplot(1,2,2);
plot(v,L_kT.','LineWidth',5);
title("$L_{kT}$ al variare di $v$");
xlabel("$v$");ylabel("$L_{kT}$ [m]");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  u_c

u_c = zeros(szu,szv,length(z));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            u_c(ii,jj,kk) = u_g(ii) + uc2_ug(ii,jj) .* ( sinh((z(kk) - z_1) ./ lambda_c) ./ sinh((z_2 - z_1) ./ lambda_c) )  + ...
                uc1_ug(ii,jj) .* ( sinh((z_2 - z(kk)) ./ lambda_c) ./ sinh((z_2 - z_1) ./ lambda_c) );
        end
    end
end

u_c_pl_z = zeros(szu,szv,length(z));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            u_c_pl_z(kk,ii,jj) = u_g(ii) + uc2_ug(ii,jj) .* ( (sinh(z(kk) ./ lambda_c)) ./ (sinh(L_g ./ lambda_c)) )  + ...
                uc1_ug(ii,jj) .* ( (sinh((L_g - z(kk)) ./ lambda_c)) ./ (sinh(L_g ./ lambda_c)) );
        end
    end
end

figure;
u_c_pl_z_pl = u_c_pl_z(:,:,3);
plot(z,u_c_pl_z_pl,'LineWidth',5);
title("$u_{c}$ al variare di $z$, con $v$ = 20");
xlabel("$z [m]$");ylabel("$u_{c}$");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  E_mn , E_cmn , E_cmn(zm)

domain = 0:0.00001:12;
bess = zeros(5,1200001);
for i = 0:2
    bess(i+1,:) = besselj(i,domain);
end

zero_0n = find(diff(sign(bess(1,:)))) * 10^-5;  
mu_0n = zero_0n(1:3);

zero_1n = find(diff(sign(bess(2,:)))) * 10^-5;  
mu_1n = zero_1n(2:4);

zero_2n = find(diff(sign(bess(3,:)))) * 10^-5;  
mu_2n = zero_2n(2:4);

clear bess; clear domain;

E_0n = ( h_bar^2 .* mu_0n.^2 ) ./ ( 2 * m_r * R.^2); 
E_0n_eV = E_0n ./ q;
E_1n = ( h_bar^2 .* mu_1n.^2 ) ./ ( 2 * m_r * R.^2);
E_1n_eV = E_1n ./ q;
E_2n = ( h_bar^2 .* mu_2n.^2 ) ./ ( 2 * m_r * R.^2); 
E_2n_eV = E_2n ./ q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_c0n = zeros(szu,szv,length(z),length(E_0n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for ll = 1:length(E_0n)
                E_c0n(ii,jj,kk,ll) = E_g/2 + E_0n(ll) - k_B .* T .* u_c(ii,jj,kk); 
            end
        end
    end
end

E_c0n_zm = zeros(szu,szv,length(E_2n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for ll = 1:length(E_0n)
            E_c0n_zm(ii,jj,ll) = E_g/2 + E_0n(ll) - k_B .* T .* u_cm(ii,jj);
        end
    end
end

E_c0n_zm_eV = E_c0n_zm ./ q;

figure;
E_c0n_zm_plot = E_c0n_zm_eV(:,:,1);
subplot(1,2,1);
plot(u_g,E_c0n_zm_plot,'linewidth',5);title("$E_{c01}(z_m)$ al variare di $u_g$");
xlabel("$u_g$");ylabel("$E_{c01}(z_m)$ [eV]");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
subplot(1,2,2);
plot(v,E_c0n_zm_plot.','linewidth',5);title("$E_{c01}(z_m)$ al variare di $v$");
xlabel("$v$");ylabel("$E_{c01}(z_m)$ [eV]");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on

figure;
E_c0n_zm_plot = E_c0n_zm_eV(:,:,2);
subplot(1,2,1);
plot(u_g,E_c0n_zm_plot,'linewidth',5);title("$E_{c02}(z_m)$ al variare di $u_g$");
xlabel("$u_g$");ylabel("$E_{c02}(z_m)$ [eV]");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
subplot(1,2,2);
plot(v,E_c0n_zm_plot.','linewidth',5);title("$E_{c02}(z_m)$ al variare di $v$");
xlabel("$v$");ylabel("$E_{c02}(z_m)$ [eV]");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_c1n = zeros(szu,szv,length(z),length(E_1n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for ll = 1:length(E_1n)
                E_c1n(ii,jj,kk,ll) = E_g/2 + E_1n(ll) - k_B .* T .* u_c(ii,jj,kk); 
            end
        end
    end
end

E_c1n_zm = zeros(szu,szv,length(E_2n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for ll = 1:length(E_1n)
            E_c1n_zm(ii,jj,ll) = E_g/2 + E_1n(ll) - k_B .* T .* u_cm(ii,jj);
        end
    end
end

E_c1n_zm_eV = E_c1n_zm ./ q;

figure;
E_c1n_zm_plot = E_c1n_zm_eV(:,:,1);
subplot(1,2,1);
plot(u_g,E_c1n_zm_plot,'linewidth',5);title("$E_{c11}(z_m)$ al variare di $u_g$");
xlabel("$u_g$");ylabel("$E_{c11}(z_m)$ [eV]");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
subplot(1,2,2);
plot(v,E_c1n_zm_plot.','linewidth',5);title("$E_{c11}(z_m)$ al variare di $v$");
xlabel("$v$");ylabel("$E_{c11}(z_m)$ [eV]");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_c2n = zeros(szu,szv,length(z),length(E_2n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for ll = 1:length(E_2n)
                E_c2n(ii,jj,kk,ll) = E_g/2 + E_2n(ll) - k_B .* T .* u_c(ii,jj,kk); 
            end
        end
    end
end

E_c2n_zm = zeros(szu,szv,length(E_2n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for ll = 1:length(E_2n)
            E_c2n_zm(ii,jj,ll) = E_g/2 + E_2n(ll) - k_B .* T .* u_cm(ii,jj);
        end
    end
end

E_c2n_zm_eV = E_c2n_zm ./ q;

figure;
E_c2n_zm_plot = E_c2n_zm_eV(:,:,1);
subplot(1,2,1);
plot(u_g,E_c2n_zm_plot,'linewidth',5);title("$E_{c21}(z_m)$ al variare di $u_g$");
xlabel("$u_g$");ylabel("$E_{c21}(z_m)$ [eV]");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
subplot(1,2,2);
plot(v,E_c2n_zm_plot.','linewidth',5);title("$E_{c21}(z_m)$ al variare di $v$");
xlabel("$v$");ylabel("$E_{c21}(z_m)$ [eV]");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% r(E)

gamma_f_z0 = zeros(szu,szv,length(epsilon_0));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for hh = 1:length(epsilon_0)
            gamma_f_z0(ii,jj,hh) = ( L_kT(ii,jj) ./ (sqrt(epsilon_0(hh)) .* lambda_0) ) .* ...
                (  atan( (L_g - z_m(ii,jj)) ./ (sqrt(epsilon_0(hh)) .* L_kT(ii,jj)) )   -  ...
                atan( ( -z_m(ii,jj)) ./ (sqrt(epsilon_0(hh)) .* L_kT(ii,jj)) )  );
        end
    end
end

r_E = gamma_f_z0 ./ (1 + gamma_f_z0);

figure;
r_E_plot = r_E(:,:,2);
subplot(1,2,1);
plot(u_g,r_E_plot,'linewidth',5);title("$r(E)$ al variare di $u_g$, $E_0 = 25 \,\, meV$");
xlabel("$u_g$");ylabel("$r(E)$");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
subplot(1,2,2);
plot(v,r_E_plot.','linewidth',5);title("$r(E)$ al variare di $v$, $E_0 = 25 \,\, meV$");
xlabel("$v$");ylabel("$r(E)$");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on

figure;
r_E_plot = r_E(:,:,9);
subplot(1,2,1);
plot(u_g,r_E_plot,'linewidth',5);title("$r(E)$ al variare di $u_g$, $E_0 = 200 \,\, meV$");
xlabel("$u_g$");ylabel("$r(E)$");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
subplot(1,2,2);
plot(v,r_E_plot.','linewidth',5);title("$r(E)$ al variare di $v$, $E_0 = 200 \,\, meV$");
xlabel("$v$");ylabel("$r(E)$");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% r(E)_mean

mu_FS = E_FS ./ (k_B*T);
mu_FD = (E_FS - E_FD) ./ (k_B*T);
% mu_FD = mu_FS - v;

eta_c0n_zm = E_c0n_zm ./ (k_B * T);
% eta_c0n_zm = zeros(szu,szv,length(E_0n));
% for ii = 1:length(u_g)
%     for jj = 1:length(v)
%         for ll = 1:length(E_0n)
%             eta_c0n_zm(ii,jj,ll) = -u_cm(ii,jj) + E_0n(ll) ./ (k_B * T);
%         end
%     end
% end
E_c0n_z0 = E_c0n(:,:,1,:);  eta_c0n_z0 = E_c0n_z0 ./ (k_B * T);

eta_c1n_zm = E_c1n_zm ./ (k_B * T);
% eta_c1n_zm = zeros(szu,szv,length(E_1n));
% for ii = 1:length(u_g)
%     for jj = 1:length(v)
%         for ll = 1:length(E_1n)
%             eta_c1n_zm(ii,jj,ll) = -u_cm(ii,jj) + E_1n(ll) ./ (k_B * T);
%         end
%     end
% end
E_c1n_z0 = E_c1n(:,:,1,:);  eta_c1n_z0 = E_c1n_z0 ./ (k_B * T);

eta_c2n_zm = E_c2n_zm ./ (k_B * T);
% eta_c2n_zm = zeros(szu,szv,length(E_2n));
% for ii = 1:length(u_g)
%     for jj = 1:length(v)
%         for ll = 1:length(E_2n)
%             eta_c2n_zm(ii,jj,ll) = -u_cm(ii,jj) + E_2n(ll) ./ (k_B * T);
%         end
%     end
% end

denom_01 = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        denom_01(ii,jj) = exp(mu_FS - eta_c0n_zm(ii,jj,1));
    end
end
denom_02 = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        denom_02(ii,jj) = exp(mu_FS - eta_c0n_zm(ii,jj,2));
    end
end
denom_11 = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        denom_11(ii,jj) = exp(mu_FS - eta_c1n_zm(ii,jj,1));
    end
end
denom_21 = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        denom_21(ii,jj) = exp(mu_FS - eta_c2n_zm(ii,jj,1));
    end
end

num_01 = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for hh = 2:length(epsilon_0)
            num_01(ii,jj) = num_01(ii,jj) + (  ( r_E(ii,jj,hh) ./ exp(epsilon_0(hh) + eta_c0n_zm(ii,jj,1)) ) .* 0.2  ); % - eta_c0n_z0
        end
    end
end
r_E_mean01 = num_01 ./ denom_01;

num_02 = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for hh = 2:length(epsilon_0)
            num_02(ii,jj) = num_02(ii,jj) + (  ( r_E(ii,jj,hh) ./ exp(epsilon_0(hh) + eta_c0n_zm(ii,jj,2)) ) .* 0.2  ); % - eta_c0n_z0(ii,jj,1,2)            
        end
    end
end
r_E_mean02 = num_02 ./ denom_02;

num_11 = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for hh = 2:length(epsilon_0)
            num_11(ii,jj) = num_11(ii,jj) + (  ( r_E(ii,jj,hh) ./ exp(epsilon_0(hh) + eta_c1n_zm(ii,jj,1)) ) .* 0.2  ); % - eta_c0n_z0            
        end
    end
end
r_E_mean11 = num_11 ./ denom_11;

num_21 = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for hh = 2:length(epsilon_0)
            num_21(ii,jj) = num_21(ii,jj) + (  ( r_E(ii,jj,hh) ./ exp(epsilon_0(hh) + eta_c2n_zm(ii,jj,1)) ) .* 0.2  ); % - eta_c0n_z0            
        end
    end
end
r_E_mean21 = num_21 ./ denom_21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Imn

I_01p = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
            I_01p(ii,jj) = (2*q*k_B*T / h) .* (1 - r_E_mean01(ii,jj)) .* exp(mu_FS - eta_c0n_zm(ii,jj,1));
%             I_01p(ii,jj) = (2*q*k_B*T / h)  .* exp(mu_FS - eta_c0n_zm(ii,jj,1));

    end
end
I_01n = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
            I_01n(ii,jj) = (2*q*k_B*T / h) .* (1 - r_E_mean01(ii,jj)) .* exp(mu_FD(jj) - eta_c0n_zm(ii,jj,1));
%             I_01n(ii,jj) = (2*q*k_B*T / h)  .* exp(mu_FD(jj) - eta_c0n_zm(ii,jj,1));

    end
end

I_01 = I_01p - I_01n;

% I_01 = zeros(szu,szv);
% for ii = 1:length(u_g)
%     for jj = 1:length(v)
%         I_01(ii,jj) = (2.*q.*k_B.*T ./ h) .* (1 - r_E_mean01(ii,jj)) .* exp( -eta_c0n_zm(ii,jj,1) ) .* (1 - exp(mu_FD(jj))); 
%     end
% end


figure;
subplot(1,2,1);
semilogy(u_g,I_01,'linewidth',5);title("$I_{01}$ al variare di $u_g$");
xlabel("$u_g$");ylabel("log($I_{01}$) $[A]$");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on;
subplot(1,2,2);
plot(v,I_01.','linewidth',5);title("$I_{01}$ al variare di $v$");
xlabel("$v$");ylabel("$I_{01}$ $[A]$");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on;

%%%%%%%%%%%%%%%%%

I_02p = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
            I_02p(ii,jj) = (2*q*k_B*T / h) .* (1 - r_E_mean02(ii,jj)) .* exp(mu_FS - eta_c0n_zm(ii,jj,2));
%             I_02p(ii,jj) = (2*q*k_B*T / h)  .* exp(mu_FS - eta_c0n_zm(ii,jj,2));

    end
end
I_02n = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
            I_02n(ii,jj) = (2*q*k_B*T / h) .* (1 - r_E_mean02(ii,jj)) .* exp(mu_FD(jj) - eta_c0n_zm(ii,jj,2));
%                 I_02n(ii,jj) = (2*q*k_B*T / h) .* exp(mu_FD(jj) - eta_c0n_zm(ii,jj,2));

    end
end

I_02 = I_02p - I_02n;

figure;
subplot(1,2,1);
semilogy(u_g,I_02,'linewidth',5);title("$I_{02}$ al variare di $u_g$");
xlabel("$u_g$");ylabel("log($I_{02}$) $[A]$");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on;
subplot(1,2,2);
plot(v,I_02.','linewidth',5);title("$I_{02}$ al variare di $v$");
xlabel("$v$");ylabel("$I_{02}$ $[A]$");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on;

%%%%%%%%%%%%%%%%%%%%%

I_11p = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
            I_11p(ii,jj) = (2*q*k_B*T / h) .* (1 - r_E_mean11(ii,jj)) .* exp(mu_FS - eta_c1n_zm(ii,jj,1));
%               I_11p(ii,jj) = (2*q*k_B*T / h) .* exp(mu_FS - eta_c1n_zm(ii,jj,1));

    end
end
I_11n = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
            I_11n(ii,jj) = (2*q*k_B*T / h) .* (1 - r_E_mean11(ii,jj)) .* exp(mu_FD(jj) - eta_c1n_zm(ii,jj,1));
%                         I_11n(ii,jj) = (2*q*k_B*T / h)  .* exp(mu_FD(jj) - eta_c1n_zm(ii,jj,1));

    end
end

I_11 = I_11p - I_11n;

figure;
subplot(1,2,1);
semilogy(u_g,I_11,'linewidth',5);title("$I_{11}$ al variare di $u_g$");
xlabel("$u_g$");ylabel("log($I_{11}$) $[A]$");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on;
subplot(1,2,2);
plot(v,I_11.','linewidth',5);title("$I_{11}$ al variare di $v$");
xlabel("$v$");ylabel("$I_{11}$ $[A]$");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_21p = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
            I_21p(ii,jj) = (2*q*k_B*T / h) .* (1 - r_E_mean21(ii,jj)) .* exp(mu_FS - eta_c2n_zm(ii,jj,1));
%                         I_21p(ii,jj) = (2*q*k_B*T / h)  .* exp(mu_FS - eta_c2n_zm(ii,jj,1));

    end
end
I_21n = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
            I_21n(ii,jj) = (2*q*k_B*T / h) .* (1 - r_E_mean21(ii,jj)) .* exp(mu_FD(jj) - eta_c2n_zm(ii,jj,1));
%                         I_21n(ii,jj) = (2*q*k_B*T / h)  .* exp(mu_FD(jj) - eta_c2n_zm(ii,jj,1));

    end
end

I_21 = I_21p - I_21n;

figure;
subplot(1,2,1);
semilogy(u_g,I_21,'linewidth',5);title("$I_{21}$ al variare di $u_g$");
xlabel("$u_g$");ylabel("log($I_{21}$) $[A]$");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on;
subplot(1,2,2);
plot(v,I_21.','linewidth',5);title("$I_{21}$ al variare di $v$");
xlabel("$v$");ylabel("$I_{21}$ $[A]$");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% I_ds

I_ds = I_01 + I_02 + I_11 + I_21;

figure;
subplot(1,2,1);
semilogy(u_g,I_ds,'linewidth',5);title("$I_{DS}$ al variare di $u_g$");
xlabel("$u_g$");ylabel("log($I_{DS}) \,\, [A]$");
legend("$v = 0$","$v = 10$","$v = 20$","$v = 30$","$v = 40$",'location','best');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
subplot(1,2,2);
plot(v,I_ds.','linewidth',5);title("$I_{DS}$ al variare di $v$");
xlabel("$v$");ylabel("$I_{DS} \,\, [A]$");
legend("$u_g = 0$","$u_g = 3.75$","$u_g = 7.5$","$u_g = 11.25$","$u_g = 15$",'location','northwest');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% f+ , f-

epsilon_0f = [5 10 15 20 25 30];

gamma_f = zeros(szu,szv,length(z),length(epsilon_0f));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for hh = 1:length(epsilon_0f)
                gamma_f(kk,hh,ii,jj) = ( L_kT(ii,jj) ./ (sqrt(epsilon_0f(hh)) .* lambda_0) ) .* ...
                    ( atan( (L_g - z_m(ii,jj)) ./ (sqrt(epsilon_0f(hh)) .* L_kT(ii,jj)) ) - ...
                    atan((z(kk) - z_m(ii,jj)) ./ (sqrt(epsilon_0f(hh)) .* L_kT(ii,jj))) );
            end
        end
    end
end

one_over_Delta_f01 = zeros(szu,szv,length(epsilon_0f));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for hh = 1:length(epsilon_0f)
            one_over_Delta_f01(ii,jj,hh) = (1 + exp(-(epsilon_0f(hh) + eta_c0n_zm(ii,jj,1) - eta_c0n_z0(ii,jj,1,1) - mu_FS))) ...
                .* (1 + gamma_f_z0(ii,jj,hh)) ;
        end
    end
end

f_p = zeros(szu,szv,length(z),length(epsilon_0f));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for hh = 1:length(epsilon_0f)
                f_p(kk,hh,ii,jj) = gamma_f(kk,hh,ii,jj) ./ one_over_Delta_f01(ii,jj,hh);
            end
        end
    end
end

f_n = f_p .* (gamma_f ./ (1 + gamma_f));

figure;
f_p_plot = f_p(:,:,4,4);
plot(z,f_p_plot,'linewidth',5);title("$f^+(E,z)$ lungo $z$, con $u_g = 11.25$ e $v = 30$, nella banda 01");
xlabel("$z$ [m]");ylabel("$f^+(E,z)$");
legend("$\epsilon_0 = 5$","$\epsilon_0 = 10$","$\epsilon_0 = 15$","$\epsilon_0 = 20$","$\epsilon_0 = 25$","$\epsilon_0 = 30$",'location','northeast');
ax = gca;
ax.FontSize = 24;
ay = gca;
ay.FontSize = 24;
grid on

figure;
f_n_plot = f_n(:,:,4,4);
plot(z,f_n_plot,'linewidth',5);title("$f^-(E,z)$ lungo $z$, con $u_g = 11.25$ e $v = 30$, nella banda 01");
xlabel("$z$ [m]");ylabel("$f^-(E,z)$");
legend("$\epsilon_0 = 5$","$\epsilon_0 = 10$","$\epsilon_0 = 15$","$\epsilon_0 = 20$","$\epsilon_0 = 25$","$\epsilon_0 = 30$",'location','northeast');
ax = gca;
ax.FontSize = 24;
ay = gca;
ay.FontSize = 24;
grid on
