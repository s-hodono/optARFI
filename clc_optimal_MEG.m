% Optimization of Encoding Gradients for MR-ARFI
% Author: Shota Hodono
% Date: Sep 2026
% fixed bug: Oct 2026

clear all
clc
% =========================================================================
%% your scanner
% =========================================================================

Gmax_yoursystem  = 200; %[mT/m]
Smax_yoursystem  = 200; %[T/m/s]


% =========================================================================
%% constants
% =========================================================================
SNR   = 20;                  % [a.u.]
gamma = 2*pi*42.577*10^6;    % [rad Hz/T]
d     = 1e-6;                % displacement [m]
ADC   = 0.7e-3 * 1e-6;       % gray matter diffusion coefficient, [m^2/s]
T2    = 80e-3;               % 80ms; T2 for gray matter at 3T [s]

Gmax = Gmax_yoursystem; %[mT/m]
Gmin = 20; %[mT/m]
S    = Smax_yoursystem; %[T/m/s]

c1   = 4e-3; % from center of excitation to start of gradient
c2   = 6e-3; % time between end of 1st lobe and start of 2nd lobe
c3   = 4e-3; % time between end of 2nd lobe to TE


% ===========================================================================
SNR_base = 5.9723; % normalization factor from following parameters
% 
% ===========================================================================
% SNR   = 20;                  % [a.u.]
% gamma = 2*pi*42.577*10^6;    % [rad Hz/T]
% d     = 1e-6;                % displacement [m]
% ADC   = 0.7e-3 * 1e-6;       % 0.7: gray matter diffusion coefficient, [m^2/s]
% T2    = 80e-3;               % T2 for gray matter at 3T [s]
% c1   = 4e-3; % from center of excitation to start of gradient
% c2   = 6e-3; % time between end of 1st lobe and start of 2nd lobe
% c3   = 4e-3; % time between end of 2nd lobe to TE
% ===========================================================================

%%
G_vec = (Gmin:0.5:Gmax)*1e-3;
t_vec = (4:.1:32)*1e-3; % total gradient duration (small delta + ramp)

% =========================================================================
%% biipolar
% =========================================================================
   
A0         = sqrt(2)*SNR*d;
for gg = 1:size(G_vec,2)
for tt = 1:size(t_vec,2)    
    G0         = G_vec(gg);
    r          = G0/S;    
    t_MEG      = t_vec(tt);
    if t_MEG/2 < 2*r 
        % triangle: G0 is not reached, slew-rate limited
        shape   = 'triangle';
        G_eff   = S*t_MEG/4;      % actual peak
        r_eff   = t_MEG/4;        % ramp = half of one lobe
        delta_b = r_eff;          % small delta = flat top + r, flat top = 0        
    else
        % trapezoid
        shape   = 'trapezoid';
        G_eff   = G0;
        r_eff   = r;
        delta_b = t_MEG/2 - r;
    end
        
    area_lobe = 2*G_eff*delta_b;

    b_val      = gamma^2*G_eff^2  * (4/3*delta_b^3 + 2*delta_b^2*r - 1/3*delta_b*r^2 + 1/15*r^3);
    TE         = c1 + c2 + c3 + 2*t_MEG;
    A2         = exp(-TE/T2);
    
    Sbp(gg,tt) = A0 * 2*gamma*area_lobe*exp(-b_val*ADC)*A2;
    
end
end

% =========================================================================
%% uniipolar
% =========================================================================

for gg = 1:size(G_vec,2)
for tt = 1:size(t_vec,2)
    G0         = G_vec(gg);
    r          = G0/S;    
    t_MEG      = t_vec(tt);
    if t_MEG < 2*r 
        r_eff   = t_MEG/2;
        G_eff   = S*r_eff;      % actual peak
        delta_u = r_eff;        % flat top = 0
    else
        % trapezoid
        r_eff   = r;
        G_eff   = G0;
        delta_u = t_MEG - r;    % flat top + r
    end
    delta_u    = t_MEG - r;
    Delta      = t_MEG + c2;
    area_lobe  = G_eff*delta_u;
    
    b_val      = gamma^2*G_eff^2 * delta_u^2* (Delta - 1/3*delta_u - 1/6*r^2/delta_u + 1/30*r^3/delta_u^2);
    TE         = c1 + c2 + c3 + 2*t_MEG;
    A2         = exp(-TE/T2);


    Suni(gg,tt) = A0 * 2*gamma*area_lobe*exp(-b_val*ADC)*A2;

end
end


% =========================================================================
%% results
% =========================================================================
disp('==================================================================')

if max(Sbp(:)) > max(Suni(:))
    disp(['use bipolar'])
    disp(['gradient shape: ', shape])
    Smax = max(Sbp(:));
else
    disp(['use unipolar'])
    Smax = max(Suni(:));
end
% Find the maximum value and its linear index
[maxVal, linearIdx] = max(Sbp(:));

% Convert linear index to row and column indices
[rowIdx, colIdx] = ind2sub(size(Sbp), linearIdx);

disp('optimal parameters')

if max(Sbp(:)) > max(Suni(:))
    disp(['G0       = ',num2str(1e3*G_vec(rowIdx)),' mT/m'])
    disp(['S        = ',num2str(Smax_yoursystem),' m/T/s'])
    disp(['t_MEG    = ',num2str(1e3*t_vec(colIdx)),' ms '])
    disp(['delta_bi = ',num2str(1e3*(t_vec(colIdx)/2 - G_vec(rowIdx)/Smax_yoursystem )),' ms '])
    disp([''])
    disp(['*** t_MEG is a total duration of two gradient lobes ***'])
else
    disp(['G0        = ',num2str(1e3*G_vec(rowIdx)),'mT/m'])
    disp(['S         = ',num2str(Smax_yoursystem),' m/T/s'])
    disp(['t_MEG     = ',num2str(1e3*t_vec(colIdx)),' ms '])
    disp(['delta_uni = ',num2str(1e3*(t_vec(colIdx) - G_vec(rowIdx)/Smax_yoursystem )),' ms '])
end
    
disp('==================================================================')

Gain = Smax / SNR_base;

disp(['relative gain vs (bipolar:G45mT/m S 200m/T/s )= ', num2str(Gain)])

% =========================================================================
%% plot
% =========================================================================
limMax = 1.2*max(Sbp(:));
figure
subplot(121)
imagesc(1e3*t_vec,(1e3*G_vec),(Sbp),[ 0 limMax])
hold on
set(gca,'YDir','normal') 
xlabel('gradient duration for two lobes [ms]','fontsize',20)
ylabel('gradient amplitude [mT/m]','fontsize',20)
title('bipolar, displacement SNR [a.u.]','fontsize',20)
grid on
colorbar
subplot(122)
imagesc(1e3*t_vec,(1e3*G_vec),(Suni),[ 0 limMax])
hold on
set(gca,'YDir','normal') 
xlabel('gradient duration for one lobe [ms]','fontsize',20)
ylabel('gradient amplitude [mT/m]','fontsize',20)
title('unipolar, displacement SNR [a.u.]','fontsize',20)
colorbar
grid on