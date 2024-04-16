% Directed Influence Model
% ========================================================================= 
% Author:          Wan-Chi Hsin
% Date:            April 16, 2024
% MATLAB Version:  R2021b
% =========================================================================
% simulation study
% single rhythm
% 4 electrodes
% 3 switching states
% indep | E1 influences E2, E4 influences E3 | E1 influences E2, E1 influences E3, E2 influences E4
% network structure is reflected through the state transition matrix (A)


%% Setting
rng(2238)
addpath(genpath(pwd))

% ===== tuning parameter =====
var_obs_noi = 5;        %var(observation noise)
aAR = 0.8;              %AR asso. w/ A

numNodes = 4;           %number of electrodes
numOsci = numNodes;     %number of oscillators
numStates = 3;          %number of switching states

fs = 100;               %sampling frequency
osciFreq = 7;           %oscillation frequency
dt = 1/fs;              %sampling interval
T = 300; TT=300;        %total time
t = (dt:dt:T);          %time axis
N = T/dt;               %number of points in data
var_proc_noi = 1;       %var(process noise)

n = numNodes;
k = numOsci;
M = numStates;

% ===== oscillation frequency =====
osc_freq = osciFreq * ones(1,k);

% ====== switch times =====
% switches at (second)
jumpat1 = 80;        
jumppoint1 = jumpat1*fs;  

jumpat2 = 200;  
jumppoint2 = jumpat2*fs;  

% true state
t1 = 0+1/fs:1/fs:jumpat1; 
t2 = jumpat1+1/fs:1/fs:jumpat2;  
t3 = jumpat2+1/fs:1/fs:TT;  
state1 = ones(size(t1));
state2 = 2*ones(size(t2));
state3 = 3*ones(size(t3));
time = [t1 t2 t3];
true_state = [state1 state2 state3];

% ===== noise =====                                                                       
% ** process noise: Q matrix **
noise = zeros(2*k,N,M);
sigma = 1;
for i=1:M
    for j=1:N
        noise(:,j,i) = sigma * mvnrnd(zeros(2*k,1), eye(2*k));
    end
end

% ** observation noise: R matrix **
obs_noi = zeros(n, N);
for i=1:n
    obs_noi(i,:) = sqrt(var_obs_noi) .* randn(1,N);
end

% ===== initial phase offsets =====
initial_phase = 2*pi*rand(1, k);

% ===== B matrix =====
% make B1 = B2 = B3
B = zeros(M,n,2*k);
for i=1:k
     B(1,i,2*i-1:2*i) = 1/sqrt(2);
end
B(2,:,:) = B(1,:,:);
B(3,:,:) = B(1,:,:);
disp(squeeze(B(1,:,:)));

% **************** A matrix ****************
alpha = zeros(n,n,M);
theta = zeros(n,n,M);
beta = zeros(n,n,M);

% linkage - preferred phase (theta) & amplitude (alpha) 
% desired amplitude scaling (beta)
% state2
theta(2,1,2)=pi/2; alpha(2,1,2)=.4; beta(2,1,2)=1;
theta(3,4,2)=-pi/2; alpha(3,4,2)=.5; beta(3,4,2)=1;

% state3
theta(2,1,3)=0; alpha(2,1,3)=.2; beta(2,1,3)=1;
theta(3,1,3)=0; alpha(3,1,3)=.4; beta(3,1,3)=1;
theta(3,2,3)=0; alpha(3,2,3)=0; beta(3,2,3)=1;
theta(4,2,3)=0; alpha(4,2,3)=.4; beta(4,2,3)=1;

% row sum of alpha
row_sum = squeeze(sum(alpha,2));
sub_A = funcs2.create_A3dim(osc_freq,fs,M,alpha,theta,beta,row_sum,aAR);

% **************** generating signals ****************
x = funcs2.carrier_signal3dim(k,N,initial_phase,jumppoint1,jumppoint2,sub_A,noise);
y = funcs2.obs_signal3(x,B,n,N,jumppoint1,jumppoint2,obs_noi);



%% Coherogram - multitaper
% *** make sure you have chronux installed ***
% specify pair to look 
iIdx = (1:4);
jIdx = (1:4);

W = 1;
T = 2;
TW = T*W;                        %time-bandwidth product given 2W=2, 
ntapers = 2*TW-1;                %... 2W is the desired freq resolution
params.Fs = 1/dt;               
params.tapers = [TW,ntapers];    
params.pad = -1;       
params.err=[2 0.05];
movingwin=[1 1];                %set the moving window dimensions 

data = y;

figure;
tiledlayout(3,2);  

pairIndex = 1;  
for i = iIdx  
    for j = jIdx  
        if i<j
            % Calculate coherogram for the current pair of electrodes
            [C,Cphi,Ccs,C1,C2,t,f] = cohgramc(data(i,:)', data(j,:)', movingwin, params); 
            % Plotting the coherence matrix
            nexttile;
            plot_matrix(C, t, f, 'n');
            colorbar;
            axis xy;
            caxis([0 1]);  
            xlabel('Time [s]');
            ylabel('Frequency [Hz]');
            title(['electrode ' num2str(i) ' & ' num2str(j)]);
            
            pairIndex = pairIndex + 1; % Update pair index
        end
    end
end
sgtitle('Coherogram');


% pair: 1 & 2
[C12,Cphi,Ccs,C1,C2,t,f] = cohgramc(data(1,:)', data(2,:)', movingwin, params); 
figure
plot_matrix(C12, t, f, 'n');
colorbar;
axis xy;
%caxis([0 1]);
xlabel('Time [s]');
ylabel('Frequency [Hz]');
title('Coherogram: Electrodes 1 & 2');


%% Theoretical cross-spectrum & coherence matrix for true A

% Initialize the coherence matrix for all states
Coh = zeros(numNodes, numNodes, numStates);
S_y = zeros(numNodes, length(f), numStates);
S_y1y2 = zeros(numNodes, numNodes,length(f), numStates);
BMat = squeeze(B(1,:,:));

for st = 1:numStates
    Aj = sub_A(:, :, st);
    Q = var_proc_noi * eye(2*numOsci);
    R = var_obs_noi * eye(numNodes);

    % Convert frequency to angular frequency in radians/sample
    omega = 2*pi*f/fs; 

    Hyy_f = zeros(numNodes, numNodes, length(f));
    Hy_f = zeros(numNodes, length(f));
    Hxx_f = zeros(2*numOsci, 2*numOsci, length(f));
    Hx_f = zeros(2*numOsci, length(f));

    % Compute theoretical cross-spectrum
    for idx = 1:length(omega)
        inv_term = inv(eye(size(Aj)) - Aj * exp(-1i*omega(idx)));
        Hxx_f(:, :, idx) = (1/fs) * inv_term * Q * inv_term';
        Hyy_f(:, :, idx) = BMat * Hxx_f(:, :, idx) * BMat' + (1/fs) * R;
        Hx_f(:, idx) = diag(Hxx_f(:, :, idx));
        Hy_f(:, idx) = diag(Hyy_f(:, :, idx));
    end

    % Compute theoretical coherence
    [~, freqIndex] = min(abs(f - osciFreq));
    S_y1y2(:,:,:,st) = Hyy_f;
    S_y(:,:,st) = Hy_f; 
    coherenceMatrix = zeros(numNodes, numNodes);

    % Compute coherence for each electrode pair at the desired frequency
    for i = 1:numNodes
        for j = 1:numNodes
            % Extract the cross-spectrum and power spectra for the electrode pair
            S_ij = S_y1y2(i, j, freqIndex,st);
            P_ii = S_y(i, freqIndex, st);
            P_jj = S_y(j, freqIndex, st);
            coherenceMatrix(i, j) = abs(S_ij) / sqrt(P_ii * P_jj);
        end
    end
    Coh(:, :, st) = abs(coherenceMatrix);
end

% plot coherence matrix for each state
upp = max(Coh,[],'all');
figure
tiledlayout(1,numStates);
for i=1:numStates
    nexttile
    imagesc(Coh(:,:,i))
    yticks(1:n); xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colorbar
    caxis([0 upp])
    colormap(parula)
    title(sprintf('A_%g',i))
end
sgtitle('Coherence (True A)')


%% Theoretical coherogram for true A
coherence = zeros(numNodes, numNodes, length(f), numStates);
for state = 1:numStates
    for i = 1:numNodes
        for j = 1:numNodes
            for freq = 1:length(f)
                S_xx = S_y(i, freq, state);
                S_yy = S_y(j, freq, state);
                S_xy = S_y1y2(i, j, freq, state);
                coherence(i, j, freq, state) = abs(S_xy) / sqrt(S_xx * S_yy);
            end
        end
    end
end

% extract pair: 1 & 2
coherence_E12 = squeeze(coherence(1,2,:,:));

% Initialize the coherogram_E12 matrix
coherogram_E12_true = zeros(N, length(f));
for i = 1:N
    state = true_state(i); 
    coherogram_E12_true(i, :) = coherence_E12(:, state)';
end

% Check if all imaginary parts are essentially zero
imagParts = imag(coherogram_E12_true);
areAllImagPartsZero = all(abs(imagParts(:)) < 1e-10); 
if areAllImagPartsZero
    disp('All complex parts are zero.');
else
    disp('There are non-zero complex parts.');
end

% plot the coherogram based on the true
figure
plot_matrix(real(coherogram_E12_true), time, f, 'n');
colorbar;
axis xy;
%caxis([0 1]); 
xlabel('Time [s]');
ylabel('Frequency [Hz]');
title('Coherogram (True A)')


%% EM on A
rng(22)

n = 4;
k = 4;
M = 3;

A = zeros(2*k,2*k,M);
B = zeros(n,2*k,M);
Q = zeros(2*k,2*k,M);
R = zeros(n,n,M);

% A matrix
A(:,:,1) = funcs2.randa(k,osc_freq(1),fs);      
A(:,:,2) = funcs2.randa(k,osc_freq(1),fs);    
A(:,:,3) = funcs2.randa(k,osc_freq(1),fs); 

% B matrix
B = zeros(n,2*k,M);
for i=1:4
    B(i,2*i-1:2*i,:) = 1/sqrt(2);
end

% R matrix
R(:,:,1) = var_obs_noi*eye(n);
R(:,:,2) = var_obs_noi*eye(n);
R(:,:,3) = var_obs_noi*eye(n);

% Q initialization
Q(:,:,1) = eye(2*k);        
Q(:,:,2) = eye(2*k);   
Q(:,:,3) = eye(2*k);   

X_0 = mvnrnd(zeros(2*k,1), eye(2*k))';
C = [.9999 .0001 .0001; .0001 .9999 .0001; .0001 .0001 .9999];

iter = 20; 
[mle_A,X_RTS,SW,Q_func] = em_projA(y',1e-5,iter,A,B,Q,R,C,X_0);


%% Plot the estimated switching state
iter = size(mle_A,4);
tm = (dt:dt:TT);
sw = SW(:,:,iter+1);

figure
tiledlayout(numStates,1);
for i=1:numStates
    nexttile
    plot(tm,sw(:,i),'LineWidth',2)
    set(gca,'YTick',[0 1],'YLim',[-0.1 1.1]);
    title(sprintf('A_%g', i));
    ylabel('Probability')
end
xlabel('Time [s]')
sgtitle('Estimated moment by moment probability of switching states')



%% Theoretical cross-spectrum & coherence for MLE A

iter = size(mle_A, 4);
mleA = mle_A(:,:,:,iter);
sw = SW(:,:,iter+1);

% the most probable state at any given time
estSt = zeros(N, 1);
for i = 1:N
    [~, estSt(i)] = max(sw(i, :));
end

% reorder index for showing network modes
mapping = zeros(1,3);
for i = 1:3
    % Find the most frequent value in estSt where S_t equals i
    modeVal = mode(estSt(true_state == i));
    mapping(i) = modeVal;
end
reorder = mapping;

% Initialize the coherence matrix for all states
Coh = zeros(numNodes, numNodes, numStates);
S_y = zeros(numNodes, length(f), numStates);
S_y1y2 = zeros(numNodes, numNodes,length(f), numStates);
BMat = squeeze(B(1,:,:));

for st = 1:numStates
    Aj = mleA(:, :, st);
    Q = var_proc_noi * eye(2*numOsci);
    R = var_obs_noi * eye(numNodes);
    
    % Convert frequency to angular frequency in radians/sample
    omega = 2*pi*f/fs; 
    Hyy_f = zeros(numNodes, numNodes, length(f));
    Hy_f = zeros(numNodes, length(f));
    Hxx_f = zeros(2*numOsci, 2*numOsci, length(f));
    Hx_f = zeros(2*numOsci, length(f));

    for idx = 1:length(omega)
        inv_term = inv(eye(size(Aj)) - Aj * exp(-1i*omega(idx)));
        Hxx_f(:, :, idx) = (1/fs) * inv_term * Q * inv_term';
        Hyy_f(:, :, idx) = BMat * Hxx_f(:, :, idx) * BMat' + (1/fs) * R;
        Hx_f(:, idx) = diag(Hxx_f(:, :, idx));
        Hy_f(:, idx) = diag(Hyy_f(:, :, idx));
    end

    [~, freqIndex] = min(abs(f - osciFreq));
    S_y1y2(:,:,:,st) = Hyy_f;
    S_y(:,:,st) = Hy_f;
    coherenceMatrix = zeros(numNodes, numNodes);

    % Compute coherence for each electrode pair at the desired frequency
    for i = 1:numNodes
        for j = 1:numNodes
            % Extract the cross-spectrum and power spectra for the electrode pair
            S_ij = S_y1y2(i, j, freqIndex,st);
            P_ii = S_y(i, freqIndex, st);
            P_jj = S_y(j, freqIndex, st);
            coherenceMatrix(i, j) = abs(S_ij) / sqrt(P_ii * P_jj);
        end
    end
    Coh(:, :, st) = abs(coherenceMatrix);
end

% plot conherence matrix for each state
upp = max(Coh,[],'all');
figure
tiledlayout(1,numStates);
for i=1:numStates
    nexttile
    imagesc(Coh(:,:,reorder(i)))
    yticks(1:n); xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colorbar
    caxis([0 upp])
    colormap(parula)
    title(sprintf('A_%g',i))
end
sgtitle('Coherence (MLE A)')


%% Theoretical coherogram for MLE A

%the most probable state at any given time
estSt = zeros(N, 1);
for i = 1:N
    [~, estSt(i)] = max(sw(i, :));
end

coherence = zeros(numNodes, numNodes, length(f), numStates);
for state = 1:numStates
    for i = 1:numNodes
        for j = 1:numNodes
            for freq = 1:length(f)
                S_xx = S_y(i, freq, state);
                S_yy = S_y(j, freq, state);
                S_xy = S_y1y2(i, j, freq, state);
                coherence(i, j, freq, state) = abs(S_xy) / sqrt(S_xx * S_yy);
            end
        end
    end
end

% extract pair: 1 & 2
coherence_E12 = squeeze(coherence(1,2,:,:));

% Initialize the coherogram_E12 matrix
coherogram_E12_mleA = zeros(N, length(f));
% weighted coherence
for i = 1:N
    coherogram_E12_mleA(i, :) = coherence_E12*sw(i,:)';
end

% Check if all imaginary parts are essentially zero
imagParts = imag(coherogram_E12_mleA);
areAllImagPartsZero = all(abs(imagParts(:)) < 1e-10); 
if areAllImagPartsZero
    disp('All complex parts are zero.');
else
    disp('There are non-zero complex parts.');
end

figure
plot_matrix(real(coherogram_E12_mleA), time, f, 'n');
colorbar;
axis xy;
%caxis([0 1]); 
xlabel('Time [s]');
ylabel('Frequency [Hz]');
title('Coherogram (MLE A)')


%% Plot three coherograms together + estimated SW

% the most probable state at any given time
estSt = zeros(N, 1);
for i = 1:N
    [~, estSt(i)] = max(sw(i, :));
end
tm = (dt:dt:TT);
newSt = estSt;
for i = 1:length(reorder)
    newSt(estSt == i) = reorder(i);
end

upp = max([max(coherogram_E12_true), max(coherogram_E12_mleA)])+0.1;
fig = figure;
set(fig, 'Position', [100, 100, 800, 650]);
ax1 = subplot(4,1,1);
plot_matrix(real(coherogram_E12_true), time, f, 'n');
colorbar;
axis xy;
caxis([0 upp]); 
ylabel('True A', 'FontSize',15);xlabel(''); title('');
text('Units', 'normalized', 'Position', [-0.12, 0.95], 'String', 'B', ...
    'FontSize', 30, 'FontWeight', 'bold');

ax2 = subplot(4,1,2);
plot_matrix(C12, t, f, 'n');
colorbar;
axis xy;
caxis([0 1]); 
ylabel('Multitaper', 'FontSize',15);xlabel(''); title('');

ax3 = subplot(4,1,3);
plot_matrix(real(coherogram_E12_mleA), time, f, 'n');
colorbar;
axis xy;
caxis([0 upp]); 
xlabel('Time [s]'); title('');
ylabel('Estimated A', 'FontSize',15);
text(-0.09, 1.88, 'Frequency [Hz]', 'Units', 'normalized', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
    'FontSize',15);

ax4 = subplot(4,1,4);
plot(tm, newSt,'b', 'LineWidth',2);
hold on
plot(tm, true_state+0.1,'r', 'LineWidth',2);
legend('Estimated State', 'True State', 'Location', 'best'); 
set(gca,'YTick',[1 2 3],'YLim',[0.8 3.2]);
xlabel('Time [s]')
ylabel('State', 'FontSize',15);
text('Units', 'normalized', 'Position', [-0.12, 0.95], 'String', 'C', ...
    'FontSize', 30, 'FontWeight', 'bold');
axesArray = [ax1, ax2, ax3, ax4];
minWidth = min(cellfun(@(x) x.Position(3), num2cell(axesArray(1:3))));
for ax = axesArray
    pos = ax.Position;
    pos(3) = minWidth; 
    ax.Position = pos;
end


%% Plot MLE A - magnitude & phase

mag = zeros(n,n,M);
phase = zeros(n,n,M);
[mag(:,:,1), phase(:,:,1)] = funcs2.A_mag_phase(mle_A(:,:,reorder(1),iter));
[mag(:,:,2), phase(:,:,2)] = funcs2.A_mag_phase(mle_A(:,:,reorder(2),iter));
[mag(:,:,3), phase(:,:,3)] = funcs2.A_mag_phase(mle_A(:,:,reorder(3),iter));

% magnitude
upp = max(mag,[],'all');

figure
tiledlayout(1,numStates);
for i=1:numStates
    nexttile
    imagesc(mag(:,:,i))
    yticks(1:n); xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colorbar
    caxis([0 upp])
    colormap(parula)
    title(sprintf('A_%g',i))
end
sgtitle('Estimated A - Magnitude')

% phase
figure
tiledlayout(1,numStates);
for i=1:numStates
    nexttile
    funcs2.imagesc2(phase(:,:,i))
    axis on
    yticks(1:n);xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colorbar
    colormap(hsv)
    caxis([-pi pi])
    title(sprintf('A_%g',i))
end
sgtitle('Estimated A - Phase')


%% Plot true A - magnitude & phase

[mag_true(:,:,1), phase_true(:,:,1)] = funcs2.A_mag_phase(sub_A(:,:,1));
[mag_true(:,:,2), phase_true(:,:,2)] = funcs2.A_mag_phase(sub_A(:,:,2));
[mag_true(:,:,3), phase_true(:,:,3)] = funcs2.A_mag_phase(sub_A(:,:,3));

% magnitude
upp = max(mag_true,[],'all');
figure
tiledlayout(1,numStates);
for i=1:numStates
    nexttile
    imagesc(mag_true(:,:,i))
    yticks(1:n); xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colorbar
    caxis([0 upp])
    colormap(parula)
    title(sprintf('A_%g',i))
end
sgtitle('True A - Magnitude')

% phase
figure
tiledlayout(1,numStates);
for i=1:numStates
    nexttile
    funcs2.imagesc2(phase_true(:,:,i))
    axis on
    yticks(1:n);xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colorbar
    colormap(hsv)
    caxis([-pi pi])
    title(sprintf('A_%g',i))
end
sgtitle('True A - Phase')


%% Plot mag & phase for true A and MLE A together
upp = max(mag_true,[],'all');
fig = figure;
set(fig, 'Position', [100, 100, 520, 650]);
tiledlayout(4,numStates);
for i=1:numStates
    nexttile
    imagesc(mag_true(:,:,i))
    yticks(1:n); xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colorbar
    caxis([0 upp])
    colormap(gca, 'parula')
    title(sprintf('State %g',i), 'FontSize', 15)
    if i==1
        ylabel('True A', 'FontSize', 15)
    end
    if i~=3
        colorbar('off')
    end
end
text('Units', 'normalized', 'Position', [-3.25, 1.05], 'String', 'D', ...
    'FontSize', 30, 'FontWeight', 'bold');
upp = max(mag,[],'all');
for i=1:numStates
    nexttile
    imagesc(mag(:,:,i))
    yticks(1:n); xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colorbar
    caxis([0 upp])
    colormap(gca, 'parula')
    if i==1
        ylabel('Estimated A', 'FontSize', 15)
    end
    if i~=3
        colorbar('off')
    end
end
text(-3.05, 1.15, 'Magnitude', 'Units', 'normalized', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'FontSize', 15,...
    'VerticalAlignment', 'bottom');
for i=1:numStates
    nexttile
    funcs2.imagesc2(phase_true(:,:,i))
    axis on
    yticks(1:n);xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colormap(gca, 'hsv')
    caxis([-pi pi]);
    cb = colorbar;
    cb.Ticks = [-pi 0 pi];
    cb.TickLabels = {'-π', '0', 'π'};
    if i==1
        ylabel('True A', 'FontSize', 15)
    end
    if i~=3
        colorbar('off')
    end
end
text('Units', 'normalized', 'Position', [-3.25, 0.95], 'String', 'E', ...
    'FontSize', 30, 'FontWeight', 'bold');
for i=1:numStates
    nexttile
    funcs2.imagesc2(phase(:,:,i))
    axis on
    yticks(1:n);xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colormap(gca, 'hsv')
    caxis([-pi pi]);
    cb = colorbar;
    cb.Ticks = [-pi 0 pi];
    cb.TickLabels = {'-π', '0', 'π'};
    if i==1
        ylabel('Estimated A', 'FontSize', 15)
    end
    if i~=3
        colorbar('off')
    end
end
text(-3.05, 1.15, 'Phase', 'Units', 'normalized', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'FontSize', 15,...
    'VerticalAlignment', 'bottom');


