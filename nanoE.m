%%% Quasi-Ballistic Transport in a nW-FET   %%%
%%% Tommaso Ugolini                         %%%
%%% University of Bologna                   %%%
%%% February 2023                           %%%

clear 
clc
close all


R = 3e-9;                       L_g = 20e-9;    % = z2 - z1
t_ox = 1e-9;                    N_D = 1e20; % 10^20 is in cm^-3
epsilon_s = 11.68;              epsilon_ox = 3.9;
k_B = physconst('Boltzmann');   T = 300;    % °K
q = 1.602176634e-19;            % electron charge [Coulomb]
E_g = q .* 1.12;                % 1.12 [eV] energy gap of silicon
n_i = 1.5e10 ; %9.65e9;         % [cm^-3] intrinsic carrier concentration
h = 6.62607015e-34;             h_bar = h/(2*pi);
m_0 = 9.1093837e-31;
m_t = 0.19 * m_0;               m_l = 0.98 * m_0;

z_1 = 0;    z_2 = L_g;
lambda_0 = 38e-9;
V = 0:0.1:1.2;
V_g = 0:0.1:1;
E_0_eV = 0:0.025:0.2;      E_0 = q .* E_0_eV;  % from 0 to 0.2 [eV] 
z = (z_1 + 1e-9):1e-9:z_2;

u_F = E_g / (2 * k_B * T);

E_F = E_g/2 + k_B * T * log(N_D / n_i);

v = (q .* V) ./ (k_B .* T);
u_g = (q .* V_g) ./ (k_B .* T);

epsilon_0 = E_0 ./ (k_B * T);

m_r = (m_t + 2 * m_l * m_t) ./ 3;

gamma = (2*epsilon_s/epsilon_ox) * log(1 + (t_ox / R));

lambda_c = R/2 * sqrt(1 + gamma);

lambda_D = sqrt( (epsilon_s * k_B * T) / (q.^2 * (N_D * 1e6)) );

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
    
% figure;
% subplot(1,2,1);
% plot(u_g,uc1_ug,'LineWidth',5);
% title("$(u_{c1} - u_g)$ al variare di $u_g$");
% xlabel("$u_g$");ylabel("$u_{c1} - u_g$");
% legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on
% subplot(1,2,2);
% plot(v,uc1_ug,'LineWidth',5);
% title("$(u_{c1} - u_g)$ al variare di $v$");
% xlabel("$v$");ylabel("$u_{c1} - u_g$");
% legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on

uc2_ug = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        uc2_ug(ii,jj) = lc_ld_2 .* ( sqrt( 1 + 2 .* ld_lc_2 .* (u_F + v(jj) - u_g(ii)) ) - 1 );
    end
end

% figure;
% subplot(1,2,1);
% plot(u_g,uc2_ug,'LineWidth',5);
% title("$(u_{c2} - u_g)$ al variare di $u_g$");
% xlabel("$u_g$");ylabel("$u_{c2} - u_g$");
% legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on
% subplot(1,2,2);
% plot(v,uc2_ug,'LineWidth',5);
% title("$(u_{c2} - u_g)$ al variare di $v$");
% xlabel("$v$");ylabel("$u_{c2} - u_g$");
% legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on

S = uc2_ug ./ uc1_ug;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  z_m

z_m = z_1 + (lambda_c / 2) .* log( (exp(L_g ./ lambda_c) - S) ./ (S - exp(-L_g ./ lambda_c)));

% figure;
% subplot(1,2,1);
% plot(u_g,z_m,'LineWidth',5);
% title("$z_{m}$ al variare di $u_g$");
% xlabel("$u_g$");ylabel("$z_{m}$ [m]");
% legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on
% subplot(1,2,2);
% plot(v,z_m,'LineWidth',5);
% title("$z_{m}$ al variare di $v$");
% xlabel("$v$");ylabel("$z_{m}$ [m]");
% legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% u_cm

u_cm = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        u_cm(ii,jj) = u_g(ii) + uc2_ug(ii,jj) .* (sinh(z_m(ii,jj)) ./ lambda_c) ./ (sinh(L_g)/lambda_c) + uc1_ug(ii,jj) .* (sinh(L_g - z_m(ii,jj)) ./ lambda_c) ./ (sinh(L_g)/lambda_c);
    end
end

% figure;
% subplot(1,2,1);
% plot(u_g,u_cm,'LineWidth',5);
% title("$u_{c}(z_m)$ al variare di $u_g$");
% xlabel("$u_g$");ylabel("$u_{c}(z_m)$");
% legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on
% subplot(1,2,2);
% plot(v,u_cm,'LineWidth',5);
% title("$u_{c}(z_m)$ al variare di $v$");
% xlabel("$v$");ylabel("$u_{c}(z_m)$");
% legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L_kT

L_kT_2 = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        L_kT_2(ii,jj) = lambda_c.^2  ./ ( k_B .* T .* (u_cm(ii,jj) - u_g(ii)) ); 
    end
end

L_kT = sqrt(L_kT_2);

% figure;
% subplot(1,2,1);
% plot(u_g,L_kT,'LineWidth',5);
% title("$L_{kT}$ al variare di $u_g$");
% xlabel("$u_g$");ylabel("$L_{kT}$ [nm]");
% legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on
% subplot(1,2,2);
% plot(v,L_kT,'LineWidth',5);
% title("$L_{kT}$ al variare di $v$");
% xlabel("$v$");ylabel("$L_{kT}$ [nm]");
% legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  u_c

u_c = zeros(szu,szv,length(z));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            u_c(ii,jj,kk) = u_g(ii) + uc2_ug(ii,jj) .* ( (sinh(z(kk)) ./ lambda_c) ./ (sinh(L_g) ./ lambda_c) )   +   uc1_ug(ii,jj) .* ( (sinh(L_g - z(kk)) ./ lambda_c) ./ (sinh(L_g) ./ lambda_c) );
        end
    end
end
u_c=real(u_c);

u_c_pl = zeros(szu,length(z),szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            u_c_pl(ii,kk,jj) = u_g(ii) + uc2_ug(ii,jj) .* ( (sinh(z(kk)) ./ lambda_c) ./ (sinh(L_g) ./ lambda_c) ) + uc1_ug(ii,jj) .* ( (sinh(L_g - z(kk)) ./ lambda_c) ./ (sinh(L_g) ./ lambda_c) );
        end
    end
end
u_c_pl = real(u_c_pl);

% %plot for z (v const)
% figure; pl=1;
% for jj = [1 5 9 13]
%     subplot(2,4,pl);
%     u_c_plot = u_c_pl(:,:,jj);
%     plot(u_g,u_c_plot);title("$u_{c}$ al variare di $u_g$, con $v$ = ", v(jj));
%     xlabel("$u_g$");ylabel("$u_{c}$");
%     ax = gca;
%     ax.FontSize = 12;
%     ay = gca;
%     ay.FontSize = 12;
%     grid on
%     subplot(2,4,pl+1);
%     plot(z,u_c_plot);title("$u_{c}$ al variare di $z$, con $v$ = ", v(jj));
%     xlabel("$z$");ylabel("$u_{c}$");
%     legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V");
%     ax = gca;
%     ax.FontSize = 12;
%     ay = gca;
%     ay.FontSize = 12;
%     grid on;
%     pl=pl+2;
% end
% 
% % plots for u_g and v (z const)
% figure;
% for zz = 1:2:length(z)
%     subplot(2,length(z)/2,zz);
%     u_c_plot = u_c(:,:,zz);
%     plot(u_g,u_c_plot);title("$u_{c}$ wrt $u_g$ with $z$ = ", z(zz));
%     xlabel("$u_g$");ylabel("$u_{c}$");
%     subplot(2,length(z)/2,zz+1);
%     plot(v,u_c_plot);title("$u_{c}$ wrt $v$ with $z$ = ", z(zz));
%     xlabel("$v$");ylabel("$u_{c}$");
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  E_mn , E_cmn , E_cmn(zm)

domain = 0:0.0001:12;
bess = zeros(5,120001);
for i = 0:2
    bess(i+1,:) = besselj(i,domain);
end

zero_0n = find(diff(sign(bess(1,:)))) * 10^-4;  
mu_0n = zero_0n(1:3);

zero_1n = find(diff(sign(bess(2,:)))) * 10^-4;  
mu_1n = zero_1n(2:4);

zero_2n = find(diff(sign(bess(3,:)))) * 10^-4;  
mu_2n = zero_2n(2:4);

clear bess; clear domain;

E_0n = ( h_bar^2 .* mu_0n.^2 ) ./ ( 2 * m_r * R.^2); 
% E_0n_eV = E_0n ./ q;
E_1n = ( h_bar^2 .* mu_1n.^2 ) ./ ( 2 * m_r * R.^2);
% E_1n_eV = E_1n ./ q;
E_2n = ( h_bar^2 .* mu_2n.^2 ) ./ ( 2 * m_r * R.^2); 
% E_2n_eV = E_2n ./ q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_c0n = zeros(szu,szv,length(z),length(E_0n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for ll = 1:length(E_0n)
                E_c0n(ii,jj,kk,ll) = E_g/2 + E_0n(ll)/2 - k_B .* T .* u_c(ii,jj,kk); 
            end
        end
    end
end

% figure;pl=1;
% for zz = [1 20]
%     subplot(2,2,pl);
%     E_c0n_plot = E_c0n(:,:,zz,1);
%     plot(u_g,E_c0n_plot,'linewidth',5);title("$E_{c01}$ al variare di $u_g$, con $z$ = ", z(zz));
%     xlabel("$u_g$");ylabel("$E_{c01}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     subplot(2,2,pl+1);
%     plot(v,E_c0n_plot,'linewidth',5);title("$E_{c01}$ al variare di $v$, con $z$ = ", z(zz));
%     xlabel("$v$");ylabel("$E_{c01}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     pl=pl+2;
% end
% figure;pl=1;
% for zz = [1 20]
%     subplot(2,2,pl);
%     E_c0n_plot = E_c0n(:,:,zz,2);
%     plot(u_g,E_c0n_plot,'linewidth',5);title("$E_{c02}$ al variare di $u_g$, con $z$ = ", z(zz));
%     xlabel("$u_g$");ylabel("$E_{c02}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     subplot(2,2,pl+1);
%     plot(v,E_c0n_plot,'linewidth',5);title("$E_{c02}$ al variare di $v$, con $z$ = ", z(zz));
%     xlabel("$v$");ylabel("$E_{c02}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     pl=pl+2;
% end
% figure;pl=1;
% for zz = [1 20]
%     subplot(2,2,pl);
%     E_c0n_plot = E_c0n(:,:,zz,3);
%     plot(u_g,E_c0n_plot,'linewidth',5);title("$E_{c03}$ al variare di $u_g$, con $z$ = ", z(zz));
%     xlabel("$u_g$");ylabel("$E_{c03}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     subplot(2,2,pl+1);
%     plot(v,E_c0n_plot,'linewidth',5);title("$E_{c03}$ al variare di $v$, con $z$ = ", z(zz));
%     xlabel("$v$");ylabel("$E_{c03}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     pl=pl+2;
% end

E_c0n_zm = zeros(szu,szv,length(E_2n));
for ii = 1:length(u_g)
    for jj = 1:length(v)

            for ll = 1:length(E_0n)
                E_c0n_zm(ii,jj,ll) = E_g/2 + E_0n(ll) - k_B .* T .* u_cm(ii,jj); 
            end
    end
end

% figure;
%     subplot(3,2,1);
%     E_c0n_zm_plot = E_c0n_zm(:,:,1);
%     plot(u_g,E_c0n_zm_plot);title("$E_{c01}(z_m)$ wrt $u_g$");
%     xlabel("$u_g$");ylabel("$E_{c01}(z_m)$");
%     subplot(3,2,2);
%     plot(v,E_c0n_zm_plot);title("$E_{c01}(z_m)$ wrt $v$");
%     xlabel("$v$");ylabel("$E_{c01}(z_m)$");
%     subplot(3,2,3);
%     E_c0n_zm_plot = E_c0n_zm(:,:,2);
%     plot(u_g,E_c0n_zm_plot);title("$E_{c02}(z_m)$ wrt $u_g$");
%     xlabel("$u_g$");ylabel("$E_{c02}(z_m)$");
%     subplot(3,2,4);
%     plot(v,E_c0n_zm_plot);title("$E_{c02}(z_m)$ wrt $v$");
%     xlabel("$v$");ylabel("$E_{c02}(z_m)$");
%     subplot(3,2,5);
%     E_c0n_zm_plot = E_c0n_zm(:,:,3);
%     plot(u_g,E_c0n_zm_plot);title("$E_{c03}(z_m)$ wrt $u_g$");
%     xlabel("$u_g$");ylabel("$E_{c03}(z_m)$");
%     subplot(3,2,6);
%     plot(v,E_c0n_zm_plot);title("$E_{c03}(z_m)$ wrt $v$");
%     xlabel("$v$");ylabel("$E_{c03}(z_m)$");
    
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

% figure;pl=1;
% for zz = [1 20]
%     subplot(2,2,pl);
%     E_c1n_plot = E_c1n(:,:,zz,1);
%     plot(u_g,E_c1n_plot,'linewidth',5);title("$E_{c11}$ al variare di $u_g$, con $z$ = ", z(zz));
%     xlabel("$u_g$");ylabel("$E_{c11}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     subplot(2,2,pl+1);
%     plot(v,E_c1n_plot,'linewidth',5);title("$E_{c11}$ al variare di $v$, con $z$ = ", z(zz));
%     xlabel("$v$");ylabel("$E_{c11}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     pl=pl+2;
% end
% figure;pl=1;
% for zz = [1 20]
%     subplot(2,2,pl);
%     E_c1n_plot = E_c1n(:,:,zz,2);
%     plot(u_g,E_c1n_plot,'linewidth',5);title("$E_{c12}$ al variare di $u_g$, con $z$ = ", z(zz));
%     xlabel("$u_g$");ylabel("$E_{c12}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     subplot(2,2,pl+1);
%     plot(v,E_c1n_plot,'linewidth',5);title("$E_{c12}$ al variare di $v$, con $z$ = ", z(zz));
%     xlabel("$v$");ylabel("$E_{c12}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     pl=pl+2;
% end
% figure;pl=1;
% for zz = [1 20]
%     subplot(2,2,pl);
%     E_c1n_plot = E_c1n(:,:,zz,3);
%     plot(u_g,E_c1n_plot,'linewidth',5);title("$E_{c13}$ al variare di $u_g$, con $z$ = ", z(zz));
%     xlabel("$u_g$");ylabel("$E_{c13}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     subplot(2,2,pl+1);
%     plot(v,E_c1n_plot,'linewidth',5);title("$E_{c13}$ al variare di $v$, con $z$ = ", z(zz));
%     xlabel("$v$");ylabel("$E_{c13}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     pl=pl+2;
% end

E_c1n_zm = zeros(szu,szv,length(E_2n));
for ii = 1:length(u_g)
    for jj = 1:length(v)

            for ll = 1:length(E_1n)
                E_c1n_zm(ii,jj,ll) = E_g/2 + E_1n(ll) - k_B .* T .* u_cm(ii,jj); 
            end
    end
end

% figure;
%     subplot(3,2,1);
%     E_c1n_zm_plot = E_c1n_zm(:,:,1);
%     plot(u_g,E_c1n_zm_plot);title("$E_{c11}(z_m)$ wrt $u_g$");
%     xlabel("$u_g$");ylabel("$E_{c11}(z_m)$");
%     subplot(3,2,2);
%     plot(v,E_c1n_zm_plot);title("$E_{c11}(z_m)$ wrt $v$");
%     xlabel("$v$");ylabel("$E_{c11}(z_m)$");
%     subplot(3,2,3);
%     E_c1n_zm_plot = E_c1n_zm(:,:,2);
%     plot(u_g,E_c1n_zm_plot);title("$E_{c12}(z_m)$ wrt $u_g$");
%     xlabel("$u_g$");ylabel("$E_{c12}(z_m)$");
%     subplot(3,2,4);
%     plot(v,E_c1n_zm_plot);title("$E_{c12}(z_m)$ wrt $v$");
%     xlabel("$v$");ylabel("$E_{c12}(z_m)$");
%     subplot(3,2,5);
%     E_c1n_zm_plot = E_c1n_zm(:,:,3);
%     plot(u_g,E_c1n_zm_plot);title("$E_{c13}(z_m)$ wrt $u_g$");
%     xlabel("$u_g$");ylabel("$E_{c13}(z_m)$");
%     subplot(3,2,6);
%     plot(v,E_c1n_zm_plot);title("$E_{c13}(z_m)$ wrt $v$");
%     xlabel("$v$");ylabel("$E_{c13}(z_m)$");
  
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

% figure;pl=1;
% for zz = [1 20]
%     subplot(2,2,pl);
%     E_c2n_plot = E_c2n(:,:,zz,1);
%     plot(u_g,E_c2n_plot,'linewidth',5);title("$E_{c21}$ al variare di $u_g$, con $z$ = ", z(zz));
%     xlabel("$u_g$");ylabel("$E_{c21}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     subplot(2,2,pl+1);
%     plot(v,E_c2n_plot,'linewidth',5);title("$E_{c21}$ al variare di $v$, con $z$ = ", z(zz));
%     xlabel("$v$");ylabel("$E_{c21}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     pl=pl+2;
% end
% figure;pl=1;
% for zz = [1 20]
%     subplot(2,2,pl);
%     E_c2n_plot = E_c2n(:,:,zz,2);
%     plot(u_g,E_c2n_plot,'linewidth',5);title("$E_{c22}$ al variare di $u_g$, con $z$ = ", z(zz));
%     xlabel("$u_g$");ylabel("$E_{c22}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     subplot(2,2,pl+1);
%     plot(v,E_c2n_plot,'linewidth',5);title("$E_{c22}$ al variare di $v$, con $z$ = ", z(zz));
%     xlabel("$v$");ylabel("$E_{c22}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     pl=pl+2;
% end
% figure;pl=1;
% for zz = [1 20]
%     subplot(2,2,pl);
%     E_c2n_plot = E_c2n(:,:,zz,3);
%     plot(u_g,E_c2n_plot,'linewidth',5);title("$E_{c23}$ al variare di $u_g$, con $z$ = ", z(zz));
%     xlabel("$u_g$");ylabel("$E_{c23}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     subplot(2,2,pl+1);
%     plot(v,E_c2n_plot,'linewidth',5);title("$E_{c23}$ al variare di $v$, con $z$ = ", z(zz));
%     xlabel("$v$");ylabel("$E_{c23}$");
%         ax = gca;
%     ax.FontSize = 18;
%     ay = gca;
%     ay.FontSize = 18;
%     grid on;
%     pl=pl+2;
% end

E_c2n_zm = zeros(szu,szv,length(E_2n));
for ii = 1:length(u_g)
    for jj = 1:length(v)

            for ll = 1:length(E_2n)
                E_c2n_zm(ii,jj,ll) = E_g/2 + E_2n(ll) - k_B .* T .* u_cm(ii,jj); 
            end
    end
end

% figure;
%     subplot(3,2,1);
%     E_c2n_zm_plot = E_c2n_zm(:,:,1);
%     plot(u_g,E_c2n_zm_plot);title("$E_{c21}(z_m)$ wrt $u_g$");
%     xlabel("$u_g$");ylabel("$E_{c21}(z_m)$");
%     subplot(3,2,2);
%     plot(v,E_c2n_zm_plot);title("$E_{c21}(z_m)$ wrt $v$");
%     xlabel("$v$");ylabel("$E_{c21}(z_m)$");
%     subplot(3,2,3);
%     E_c2n_zm_plot = E_c2n_zm(:,:,2);
%     plot(u_g,E_c2n_zm_plot);title("$E_{c22}(z_m)$ wrt $u_g$");
%     xlabel("$u_g$");ylabel("$E_{c22}(z_m)$");
%     subplot(3,2,4);
%     plot(v,E_c2n_zm_plot);title("$E_{c22}(z_m)$ wrt $v$");
%     xlabel("$v$");ylabel("$E_{c22}(z_m)$");
%     subplot(3,2,5);
%     E_c2n_zm_plot = E_c2n_zm(:,:,3);
%     plot(u_g,E_c2n_zm_plot);title("$E_{c23}(z_m)$ wrt $u_g$");
%     xlabel("$u_g$");ylabel("$E_{c23}(z_m)$");
%     subplot(3,2,6);
%     plot(v,E_c2n_zm_plot);title("$E_{c23}(z_m)$ wrt $v$");
%     xlabel("$v$");ylabel("$E_{c23}(z_m)$");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma_f_z0 = zeros(szu,szv,length(epsilon_0));
for ii = 1:length(u_g)
    for jj = 1:length(v)
            for hh = 1:length(epsilon_0)
                gamma_f_z0(ii,jj,hh) = ( L_kT(ii,jj) ./ (sqrt(epsilon_0(hh)) .* lambda_0) ) .* (  atan( (L_g - z_m(ii,jj)) ./ (sqrt(epsilon_0(hh)) .* L_kT(ii,jj)) )   -   atan( ( -z_m(ii,jj)) ./ (sqrt(epsilon_0(hh)) .* L_kT(ii,jj)) )  );
            end
    end
end

r_E = gamma_f_z0 ./ (1 + gamma_f_z0);

% figure; pl=1;
% for e0e0 = [ 1 , 3 , 5 , 7 , 9 ]
%     subplot(3,4,pl);
%     r_E_plot = r_E(:,:,e0e0);
%     plot(u_g,r_E_plot);title("$r(E)$ wrt $u_g$ with $E_0$ = ", E_0_eV(e0e0));
%     xlabel("$u_g$");ylabel("$r(E)$");
%     subplot(3,4,pl + 1);
%     plot(v,r_E_plot);title("$r(E)$ wrt $v$ with $E_0$ = ", E_0_eV(e0e0));
%     xlabel("$v$");ylabel("$r(E)$");
%     pl = pl + 2;
% end

sum_r_E = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for hh = 2:length(epsilon_0)
            sum_r_E(ii,jj) = sum_r_E(ii,jj) + r_E(ii,jj,hh);
        end
    end
end

weightedsum_r_Ec0n = zeros(szu,szv,length(z),length(E_0n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for ll = 1:length(E_0n)
                for hh = 2:length(epsilon_0)
                    weightedsum_r_Ec0n(ii,jj,kk,ll) = weightedsum_r_Ec0n(ii,jj,kk,ll) + ( log(1 + exp(E_c0n(ii,jj,kk,ll) - E_F)) .* r_E(ii,jj,hh));
                end
            end
        end
    end
end
r_Ec0n_mean = zeros(szu,szv,length(z),length(E_0n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for ll = 1:length(E_0n)
                    r_Ec0n_mean(ii,jj,kk,ll) = weightedsum_r_Ec0n(ii,jj,kk,ll) ./ sum_r_E(ii,jj);
            end
        end
    end
end

weightedsum_r_Ec1n = zeros(szu,szv,length(z),length(E_1n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for ll = 1:length(E_1n)
                for hh = 2:length(epsilon_0)
                    weightedsum_r_Ec1n(ii,jj,kk,ll) = weightedsum_r_Ec1n(ii,jj,kk,ll)  +  ( log( 1 + exp(E_c1n(ii,jj,kk,ll) - E_F) ) .* r_E(ii,jj,hh) );
                end
            end
        end
    end
end
r_Ec1n_mean = zeros(szu,szv,length(z),length(E_1n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for ll = 1:length(E_1n)
                    r_Ec1n_mean(ii,jj,kk,ll) = weightedsum_r_Ec1n(ii,jj,kk,ll) ./ sum_r_E(ii,jj);
            end
        end
    end
end

weightedsum_r_Ec2n = zeros(szu,szv,length(z),length(E_2n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for ll = 1:length(E_2n)
                for hh = 2:length(epsilon_0)
                    weightedsum_r_Ec2n(ii,jj,kk,ll) = weightedsum_r_Ec2n(ii,jj,kk,ll) + ( log(1 + exp(E_c2n(ii,jj,kk,ll) - E_F)) .* r_E(ii,jj,hh));
                end
            end
        end
    end
end
r_Ec2n_mean = zeros(szu,szv,length(z),length(E_2n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for ll = 1:length(E_2n)
                    r_Ec2n_mean(ii,jj,kk,ll) = weightedsum_r_Ec2n(ii,jj,kk,ll) ./ sum_r_E(ii,jj);
            end
        end
    end
end

mu_FS = E_F ./ (k_B*T);
mu_FD = (E_F + q .* V) ./ (k_B*T);

% eta_c0n = E_c0n ./ (k_B * T);
% eta_c1n = E_c1n ./ (k_B * T);
% eta_c2n = E_c2n ./ (k_B * T);

eta_c0n_zm = E_c0n_zm ./ (k_B * T);
eta_c1n_zm = E_c1n_zm ./ (k_B * T);
eta_c2n_zm = E_c2n_zm ./ (k_B * T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Imn

I_0n = zeros(szu,szv,length(E_0n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
            for ll = 1:length(E_0n)
                    I_0n(ii,jj,ll) = (2*q*k_B*T / h) .* (1 - r_Ec0n_mean(ii,jj,1,ll)) .* log((1 + exp(mu_FS - eta_c0n_zm(ii,jj,ll)))  ./  (1 + exp(mu_FD(jj) - eta_c0n_zm(ii,jj,ll))));
            end
    end
end

% figure;
% subplot(1,2,1);
% I_0n_plot = I_0n(:,:,1);
% plot(u_g,I_0n_plot,'linewidth',5);title("$I_{01}$ al variare di $u_g$");
% xlabel("$u_g$");ylabel("$I_{01}$ $[A]$");
% legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;
% subplot(1,2,2);
% plot(v,I_0n_plot,'linewidth',5);title("$I_{01}$ al variare di $v$");
% xlabel("$v$");ylabel("$I_{01}$ $[A]$");
% legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;
% figure;
% subplot(1,2,1);
% I_0n_plot = I_0n(:,:,2);
% plot(u_g,I_0n_plot,'linewidth',5);title("$I_{02}$ al variare di $u_g$");
% xlabel("$u_g$");ylabel("$I_{02}$ $[A]$");
% legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$",'location','southwest');
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;
% subplot(1,2,2);
% plot(v,I_0n_plot,'linewidth',5);title("$I_{02}$ al variare di $v$");
% xlabel("$v$");ylabel("$I_{02}$ $[A]$");
% legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V",'location','southwest');
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;
% figure;
% subplot(1,2,1);
% I_0n_plot = I_0n(:,:,3);
% plot(u_g,I_0n_plot,'linewidth',5);title("$I_{03}$ al variare di $u_g$");
% xlabel("$u_g$");ylabel("$I_{03}$ $[A]$");
% legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;
% subplot(1,2,2);
% plot(v,I_0n_plot,'linewidth',5);title("$I_{03}$ al variare di $v$");
% xlabel("$v$");ylabel("$I_{03}$ $[A]$");
% legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V");
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;

I_1n = zeros(szu,szv,length(E_1n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
            for ll = 1:length(E_1n)
                    I_1n(ii,jj,ll) = 2*q*k_B*T / h .* (1 - r_Ec1n_mean(ii,jj,1,ll)) .* log((1 + exp(mu_FS - eta_c1n_zm(ii,jj,ll)))  ./  (1 + exp(mu_FD(jj) - eta_c1n_zm(ii,jj,ll))));
            end
    end
end

% figure;
% subplot(1,2,1);
% I_1n_plot = I_1n(:,:,1);
% plot(u_g,I_1n_plot,'linewidth',5);title("$I_{11}$ al variare di $u_g$");
% xlabel("$u_g$");ylabel("$I_{11}$ $[A]$");
% legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$",'location','southwest');
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;
% subplot(1,2,2);
% plot(v,I_1n_plot,'linewidth',5);title("$I_{11}$ al variare di $v$");
% xlabel("$v$");ylabel("$I_{11}$ $[A]$");
% legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V",'location','southwest');
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;
% figure;
% subplot(1,2,1);
% I_1n_plot = I_1n(:,:,2);
% plot(u_g,I_1n_plot,'linewidth',5);title("$I_{12}$ al variare di $u_g$");
% xlabel("$u_g$");ylabel("$I_{12}$ $[A]$");
% legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$",'location','southwest');
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;
% subplot(1,2,2);
% plot(v,I_1n_plot,'linewidth',5);title("$I_{12}$ al variare di $v$");
% xlabel("$v$");ylabel("$I_{12}$ $[A]$");
% legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V",'location','southwest');
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;

I_2n = zeros(szu,szv,length(E_2n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
            for ll = 1:length(E_2n)
                    I_2n(ii,jj,ll) = 2*q*k_B*T / h .* (1 - r_Ec2n_mean(ii,jj,1,ll)) .* log((1 + exp(mu_FS - eta_c2n_zm(ii,jj,ll)))  ./  (1 + exp(mu_FD(jj) - eta_c2n_zm(ii,jj,ll))));
            end
    end
end

% figure;
% subplot(1,2,1);
% I_2n_plot = I_2n(:,:,1);
% plot(u_g,I_2n_plot,'linewidth',5);title("$I_{21}$ al variare di $u_g$");
% xlabel("$u_g$");ylabel("$I_{21}$ $[A]$");
% legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$",'location','southwest');
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;
% subplot(1,2,2);
% plot(v,I_2n_plot,'linewidth',5);title("$I_{21}$ al variare di $v$");
% xlabel("$v$");ylabel("$I_{21}$ $[A]$");
% legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V",'location','southwest');
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;
% figure;
% subplot(1,2,1);
% I_2n_plot = I_2n(:,:,2);
% plot(u_g,I_2n_plot,'linewidth',5);title("$I_{22}$ al variare di $u_g$");
% xlabel("$u_g$");ylabel("$I_{22}$ $[A]$");
% legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$",'location','southwest');
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;
% subplot(1,2,2);
% plot(v,I_2n_plot,'linewidth',5);title("$I_{22}$ al variare di $v$");
% xlabel("$v$");ylabel("$I_{22}$ $[A]$");
% legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V",'location','southwest');
% ax = gca;
% ax.FontSize = 22;
% ay = gca;
% ay.FontSize = 22;
% grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% I_ds

I_ds = zeros(szu,szv);
for ii = 1:length(u_g)
    for jj = 1:length(v)
            for ll = 1:length(E_0n)
                    I_ds(ii,jj) = I_ds(ii,jj) + (I_0n(ii,jj,ll) + I_1n(ii,jj,ll) + I_2n(ii,jj,ll));
            end
    end
end

figure;
subplot(1,2,1);
plot(u_g,I_ds,'linewidth',5);title("$I_{DS}$ al variare di $u_g$");
xlabel("$u_g$");ylabel("$I_{DS} \,\, [A]$");
legend("$V_{DS} = 0 V$","$V_{DS} = 0.1 V$","$V_{DS} = 0.2 V$","$V_{DS} = 0.3 V$","$V_{DS} = 0.4 V$","$V_{DS} = 0.5 V$","$V_{DS} = 0.6 V$","$V_{DS} = 0.7 V$","$V_{DS} = 0.8 V$","$V_{DS} = 0.9 V$","$V_{DS} = 1 V$","$V_{DS} = 1.1 V$","$V_{DS} = 1.2 V$",'location','southwest');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on
subplot(1,2,2);
plot(v,I_ds,'linewidth',5);title("$I_{DS}$ al variare di $v$");
xlabel("$v$");ylabel("$I_{DS} \,\, [A]$");
legend("$V_g = 0$ V","$V_g = 0.1$ V","$V_g = 0.2$ V","$V_g = 0.3$ V","$V_g = 0.4$ V","$V_g = 0.5$ V","$V_g = 0.6$ V","$V_g = 0.7$ V","$V_g = 0.8$ V","$V_g = 0.9$ V","$V_g = 1$ V",'location','southwest');
ax = gca;
ax.FontSize = 22;
ay = gca;
ay.FontSize = 22;
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fermi_integral = TBFD_integral(0,'C:\Users\tommu\OneDrive\Documents\MATLAB\index.db');
% 
% f0_of_mufs_etac0n = log(1 + exp(mu_FS - eta_c0n_zm));
% f0_of_mufs_etac1n = log(1 + exp(mu_FS - eta_c1n_zm));
% f0_of_mufs_etac2n = log(1 + exp(mu_FS - eta_c2n_zm));
% 
% % r_E_mean_0 = 1 ./ ( k_B .* T .* f0_of_mufs_etac0n) .* integrate();


gamma_f = zeros(szu,szv,length(z),length(epsilon_0));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for hh = 1:length(epsilon_0)
                gamma_f(kk,hh,ii,jj) = ( L_kT(ii,jj) ./ (sqrt(epsilon_0(hh)) .* lambda_0) ) .* ( atan( (L_g - z_m(ii,jj)) ./ (sqrt(epsilon_0(hh)) .* L_kT(ii,jj)) ) - ...
                    atan((z(kk) - z_m(ii,jj)) ./ (sqrt(epsilon_0(hh)) .* L_kT(ii,jj))) );
            end
        end
    end
end

Delta_f = log(1 + exp(E_F)) - log(1 + exp(E_F)) .* gamma_f_z0;

f_p = zeros(szu,szv,length(z),length(E_2n));
for ii = 1:length(u_g)
    for jj = 1:length(v)
        for kk = 1:length(z)
            for hh = 1:length(epsilon_0)
                f_p(kk,hh,ii,jj) = Delta_f(ii,jj,hh) .* gamma_f(kk,hh,ii,jj);
            end
        end
    end
end

f_n = f_p .* (gamma_f ./ (1 + gamma_f));

% figure;
% f_p_plot = f_p(:,:,4,4);
% plot(z,f_p_plot,'linewidth',5);title("$f^+$ lungo $z$");
% xlabel("$z$");ylabel("$f^+$");
% legend("$E_0 = 0$ meV","$E_0 = 25$ meV","$E_0 = 50$ meV","$E_0 = 75$ meV","$E_0 = 100$ meV","$E_0 = 125$ meV","$E_0 = 150$ meV","$E_0 = 175$ meV","$E_0 = 200$ meV",'location','northeast');
% ax = gca;
% ax.FontSize = 24;
% ay = gca;
% ay.FontSize = 24;
% grid on
% 
% figure;
% f_n_plot = f_n(:,:,4,4);
% plot(z,f_n_plot,'linewidth',5);title("$f^-$ lungo $z$");
% xlabel("$z$");ylabel("$f^-$");
% legend("$E_0 = 0$ meV","$E_0 = 25$ meV","$E_0 = 50$ meV","$E_0 = 75$ meV","$E_0 = 100$ meV","$E_0 = 125$ meV","$E_0 = 150$ meV","$E_0 = 175$ meV","$E_0 = 200$ meV",'location','northeast');
% ax = gca;
% ax.FontSize = 24;
% ay = gca;
% ay.FontSize = 24;
% grid on





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subplot(1,2,2);
% plot(epsilon_0,f_p_plot);title("$f^+$ wrt $\epsilon_0$");
% xlabel("$\epsilon_0$");ylabel("$f^+$");


% figure;pl=1;
% for zz = 1:4:length(z)
%     subplot(2,round(length(z)/4),pl);
%     f_p_plot = f_p(:,:,zz,5);
%     plot(u_g,f_p_plot);title("$E_{c23}$ wrt $u_g$ with $z$ = ", z(zz));
%     xlabel("$u_g$");ylabel("$E_{c23}$");
%     subplot(2,round(length(z)/4),pl+1);
%     plot(v,f_p_plot);title("$E_{c23}$ wrt $v$ with $z$ = ", z(zz));
%     xlabel("$v$");ylabel("$E_{c23}$");
%     pl=pl+2;
% end
% figure;pl=1;
% for zz = 1:4:length(z)
%     subplot(2,round(length(z)/4),pl);
%     f_p_plot = f_p(:,:,zz,9);
%     plot(u_g,f_p_plot);title("$E_{c22}$ wrt $u_g$ with $z$ = ", z(zz));
%     xlabel("$u_g$");ylabel("$E_{c22}$");
%     subplot(2,round(length(z)/4),pl+1);
%     plot(v,f_p_plot);title("$E_{c22}$ wrt $v$ with $z$ = ", z(zz));
%     xlabel("$v$");ylabel("$E_{c22}$");
%     pl=pl+2;
% end

% f_p0 = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







