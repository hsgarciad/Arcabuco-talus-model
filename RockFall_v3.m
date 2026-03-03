%% Talus burial and erosion model assuming triangular shape%%
%This model builds up an inventory of atoms from which we derive a
%concentration and ratio at every time step
clear;clc
%% Input Parameters %%

tsim = 1e6;                           %Total simulation time
dt = 1;                               %Time step (years)
rp = 89;                              %Return period (years)
theta = 8;                            %Angle of repose
phi = 0.40;                           %Talus porosity
e_back = 4.83e-6;                     %Background erosion rate (m/yr)
rho_rock = 2.65;                      %Average quartz density (g/cm^3)
rho_talus = rho_rock * (1-phi);       %Effective density in the deposit
L_cliff = 26.3;                         %Average Cliff height used to scale V to Z (m)


%- Cosmogenic 10Be and 26Al constants -%
P10_0 = 18.6;                        % 10Be surface prod rate from watershed 1 (atoms/g/yr)
P26_0 = 131;                         % 26Al surface prod rate from watershed 1 (atoms/g/yr)
L10 = 5.0011e-07;                    % 10Be decay constant (yr^-1) after Chmeleff et al. 2010
L26 = 9.6673e-07;                    % 26Al decay constant (yr^-1) after Nudat 3.0 at https://www.nndc.bnl.gov/nudat3/
Lambda = 160;                        % Attenuation length (g/cm^2)

%- Pre-allocate zeros -%
H = 0.3;                            % Start with 30 cm of material - average block size in Arcabuco
V_storage = (H^3) / (2 * tand(theta)^2);
Inv10 = 0; 
Inv26 = 0;

C10_export = zeros(tsim, 1);
C26_export = zeros(tsim, 1);
% ratio_history = zeros(tsim, 1);
H_history = zeros(tsim, 1);

%% Simulation Loop
for t = 1:tsim
    % --- A. Production in Storage (Exact Triangular Solution) ---
    if V_storage > 0.001
        H = nthroot(2 * V_storage * tand(theta)^2, 3); 
                         
        
        % Mass depth calculation (converting H_meters to H_cm for Lambda)
        z_rho = H * 100 * rho_talus; 
        
        % Accounts for thinning from cliff (H) to toe (0)
        % Formula: (P0 * Lambda * W / rho) * [1 - (Lambda / z_rho) * (1 - exp(-z_rho/Lambda))]
        
        avg_P10 = 2 * (P10_0 * Lambda / z_rho) * (1 - (Lambda/z_rho) * (1 - exp(-z_rho/Lambda)));
        avg_P26 = 2 * (P26_0 * Lambda / z_rho) * (1 - (Lambda/z_rho) * (1 - exp(-z_rho/Lambda)));
        P10_total = avg_P10 * V_storage * 1e6 * rho_talus;
        P26_total = avg_P26 * V_storage * 1e6 * rho_talus;
        
        % Update Inventories (Production + Decay)
        Inv10 = Inv10 + P10_total - (Inv10 * L10);
        Inv26 = Inv26 + P26_total - (Inv26 * L26);
    end
    
    % --- B. Stochastic Rockfall (V-F Scaling & Square Slab Fix) ---
    if mod(t, rp) == 0
        % Parameters from Moos-style scaling - Deterministic
        Ast = 0.03; B = 0.5; S = 1.526; %From fracture data 
        Ft = 1 / rp;   
        % Frequency-Volume Scaling Equation
        V_event = min(((Ast * S) / Ft)^(1/B), 200);
        
        % %Quasi Stochastic
        % B_event = 0.5; Ast = 0.03; S = 10; %From fracture data 
        % Ft = 1 / rp;
        % V_char = ((Ast * S) / Ft)^(1/B_event);  % characteristic volume anchors the scale
        % sigma = 0.5;  % controls spread
        % V_event = min(V_char * exp(sigma * randn()), 200);
        
                           
        % Calculate Inheritance Depth (Square Slab Logic)
        % Z_max (m) = sqrt(Volume / Cliff_Height)
        Z_max = sqrt(V_event / L_cliff);
        z_cm = Z_max * 100;
        
        % Average Concentration of the falling block
        z_g = z_cm * rho_rock;   % g/cm²
        % Vertical erosion term — surface lowering
        eterm_vert  = rho_rock * (e_back * 100) / Lambda;

        % Horizontal erosion term — cliff retreat from rockfall flux
        retreat_cliff    = (V_event / rp) / (S * 1e4);  % m/yr
        eterm_horiz = rho_rock * (retreat_cliff * 100) / Lambda;

        % Two-surface inheritance
        C10_horiz = (P10_0/2  / (L10 + eterm_horiz)) * (Lambda/z_g) * (1 - exp(-z_g/Lambda));
        C10_vert  = (P10_0    / (L10 + eterm_vert))  * (Lambda/z_g) * (1 - exp(-z_g/Lambda));
        C10_inh   = C10_horiz + C10_vert;

        C26_horiz = (P26_0/2  / (L26 + eterm_horiz)) * (Lambda/z_g) * (1 - exp(-z_g/Lambda));
        C26_vert  = (P26_0    / (L26 + eterm_vert))  * (Lambda/z_g) * (1 - exp(-z_g/Lambda));
        C26_inh   = C26_horiz + C26_vert;

        % Update Storage Volume and Inventories
        % Note: V_event is solid rock; V_storage is bulk (including porosity)
        V_storage = V_storage + (V_event / (1 - phi));
        Inv10 = Inv10 + (V_event * 1e6 * C10_inh * rho_rock); % Mass = Vol_solid * rho_rock
        Inv26 = Inv26 + (V_event * 1e6 * C26_inh * rho_rock);
    end
    
    % --- C. Erosion and Export (Background Drawdown) ---
    if V_storage > 0
        H_old = nthroot(2 * V_storage * tand(theta)^2, 3);
        H_new = max(0, H_old - e_back); 
        
        V_new = (H_new^3) / (2 * tand(theta)^2);
        V_out = max(0, V_storage - V_new); 
        
        % Concentration of Export (Calculated before subtracting volume)
        conc10_now = Inv10 / (V_storage * 1e6 * rho_talus);
        conc26_now = Inv26 / (V_storage * 1e6 * rho_talus);
        
        % Update Atoms and Geometry
        Inv10 = max(0, Inv10 - (V_out * 1e6 * conc10_now * rho_talus));
        Inv26 = max(0, Inv26 - (V_out * 1e6 * conc26_now * rho_talus));
        
        V_storage = V_new;
        H = H_new;
        
        % Record export concentrations for analysis
        C10_export(t) = conc10_now;
        C26_export(t) = conc26_now;

        % ratio_history(t) = C26_export(t) ./ C10_export(t);
        % ratio_history(ratio_history == 0) = NaN;
    end
    
    H_history(t) = H;
end

%% 5. VISUALIZATION
figure('Position', [100 50 800 900]);

subplot(3,1,1);
plot(H_history, 'Color', [0.4 0.2 0]);
ylabel('Thickness H (m)'); title('Basin Filling & Stripping'); grid on;

subplot(3,1,2);
plot(C10_export, 'b'); hold on;
plot(C26_export, 'r--');
ylabel('Conc (at/g)'); legend('^{10}Be', '^{26}Al (norm)'); grid on;

subplot(3,1,3);
plot(C26_export ./ C10_export, 'k');
yline(P26_0/P10_0, 'r:', 'Surface Production Ratio');
xline(381e3,      'g--', 'Simple burial age', 'LineWidth', 1.5);
xline(743e3,                'm--', 'Model age', 'LineWidth', 1.5);
ylabel('Ratio ^{26}Al/^{10}Be'); xlabel('Time (Years)'); grid on;