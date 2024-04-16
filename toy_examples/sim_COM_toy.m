% Common Oscillator Model
% ========================================================================= 
% Author:          Wan-Chi Hsin
% Date:            April 16, 2024
% MATLAB Version:  R2021b
% =========================================================================
% simulation study
% single rhythm
% 4 electrodes
% 3 switching states
% indep | E1.E2 Linked & E3.E4 linked | E1.E2.E3 linked
% network structure is reflected through the observation matrix (B)


%% Setting
rng(2238)
addpath(genpath(pwd))

% ===== tuning parameter =====
var_obs_noi = 3;        %var(observation noise)
aAR = 0.85;             %AR asso. w/ A

numNodes = 4;           %number of electrodes
numOsci = 2;            %number of oscillators
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
osc_freq = [osciFreq osciFreq]; 

% ====== switch =====
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

% rho for each oscillator
rho = aAR * ones(k,M); 

% ===== noise =====                                                                     
% ** process noise: Q matrix **
sigma = sqrt(var_proc_noi);
noise = zeros(2*k,N,M);

% ** state1
e1 = zeros(2*k, jumppoint1);
for i = 1:jumppoint1
    e1(:,i) = sigma * mvnrnd(zeros(2*k,1), eye(2*k));  
end
noise(:,1:jumppoint1,1) = e1;
% ** state2
e2 = zeros(4, jumppoint2-jumppoint1);
for i = 2:jumppoint2-jumppoint1+1
    e2(:,i) = sigma * mvnrnd(zeros(2*k,1), eye(2*k));  
end
noise(:,jumppoint1+1:jumppoint2,2) = e2(:,2:jumppoint2-jumppoint1+1);
% ** state3
e3 = zeros(4, N-jumppoint2);
for i = 2:N-jumppoint2+1
    e3(:,i) = sigma * mvnrnd(zeros(2*k,1), eye(2*k)); 
end
noise(:,jumppoint2+1:N,3) = e3(:,2:N-jumppoint2+1);

% ** observation noise: R matrix **
obs_noi = zeros(n, N);
for i=1:n
    obs_noi(i,:) = sqrt(var_obs_noi) * randn(1,N);
end

% ===== initial phase offsets =====
initial_phase = 2*pi*rand(1, k);

% ===== B matrix (design here: network structure) =====
B = zeros(M,n,2*k);
B(1,1,1) = .4; B(1,2,3) = .4;

B(2,1,1) = .3; B(2,2,2) = .3; 
B(2,3,3) = .3; B(2,4,4) = .3; 

B(3,1,1) = .25; B(3,2,1) = -.25;
B(3,3,1) = .25; B(3,4,3) = .2;

disp('True B matrix');
disp(squeeze(B(1,:,:)));
disp(squeeze(B(2,:,:)));
disp(squeeze(B(3,:,:)));

% **************** generating signals ****************
sub_A = funcs2.create_A3(osc_freq,fs,rho,M);
x = funcs2.carrier_signal3(k,N,initial_phase,jumppoint1,jumppoint2,sub_A,noise);
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
title('Coherogram (multitaper)');


%%  Theoretical cross-spectrum & coherence matrix for true B

% Initialize the coherence matrix for all states
Coh = zeros(numNodes, numNodes, numStates);
S_y = zeros(numNodes, length(f), numStates);
S_y1y2 = zeros(numNodes, numNodes,length(f), numStates);

theta1 = (2*pi*osc_freq(1))*(1/fs);
mat1 = [cos(theta1) -sin(theta1); sin(theta1) cos(theta1)];
A = blkdiag(rho(1,1)*mat1, rho(2,1)*mat1);

for st = 1:numStates
    Bj = squeeze(B(st, :, :));
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
        inv_term = inv(eye(size(A)) - A * exp(-1i*omega(idx)));
        Hxx_f(:, :, idx) = (1/fs) * inv_term * Q * inv_term';
        Hyy_f(:, :, idx) = Bj * Hxx_f(:, :, idx) * Bj' + (1/fs) * R;
        Hx_f(:, idx) = diag(Hxx_f(:, :, idx));
        Hy_f(:, idx) = diag(Hyy_f(:, :, idx));
    end

    % Compute theoretical coherence
    [~, freqIndex] = min(abs(f - osciFreq));
    S_y1y2(:,:,:,st) = Hyy_f;
    S_y(:,:,st) = Hy_f; 
    coherenceMatrix = zeros(numNodes, numNodes);

    for i = 1:numNodes
        for j = 1:numNodes
            % Extract the cross-spectrum and power spectra for the electrode pair
            S_ij = S_y1y2(i, j, freqIndex, st);
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
    title(sprintf('B_%g',i))
end
sgtitle('Coherence (True B)')


%% Theoretical coherogram for true B
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
title('Coherogram (True B)')


%% EM on B
rng(60)

k = numOsci;
% initialization
A = zeros(2*k,2*k,M);
H = zeros(n,2*k,M);
Q = zeros(2*k,2*k,M);
R = zeros(n,n,M);

theta1 = (2*pi*osc_freq(1))*(1/fs);
theta2 = (2*pi*osc_freq(2))*(1/fs);
mat1 = [cos(theta1) -sin(theta1); sin(theta1) cos(theta1)];
mat2 = [cos(theta2) -sin(theta2); sin(theta2) cos(theta2)];

% state1
A(:,:,1) = blkdiag(rho(1,1)*mat1, rho(2,1)*mat2);
Q(:,:,1) = var_proc_noi * eye(2*k);
R(:,:,1) = var_obs_noi * eye(n);

% state2
A(:,:,2) = blkdiag(rho(1,2)*mat1, rho(2,2)*mat2);
Q(:,:,2) = var_proc_noi * eye(2*k);
R(:,:,2) = var_obs_noi * eye(n);

% state3
A(:,:,3) = blkdiag(rho(1,3)*mat1, rho(2,3)*mat2);
Q(:,:,3) = var_proc_noi * eye(2*k);  
R(:,:,3) = var_obs_noi * eye(n);

% B initialization
% everything is randomly weakly correlated
B0 = .05*rand(n,2*k,M);    

X_0 = mvnrnd(zeros(2*k,1), eye(2*k))';
C = [.9999 .0001 .0001; .0001 .9999 .0001; .0001 .0001 .9999];

iter = 20; 
[mle_B,X_RTS,SW,Q_func] = em_B(y',1e-5,iter,A,B0,Q,R,C,X_0);


%% Plot the estimateed switching state
iter = size(mle_B,4);
tm = (dt:dt:TT);
sw = SW(:,:,iter+1);

figure
tiledlayout(numStates,1);
for i=1:numStates
    nexttile
    plot(tm,sw(:,i),'LineWidth',2)
    set(gca,'YTick',[0 1],'YLim',[-0.1 1.1]);
    title(sprintf('B_%g', i));
    ylabel('Probability')
end
xlabel('Time [s]')
sgtitle('Estimated moment by moment probability of switching states')


%%  Theoretical cross-spectrum & coherence for MLE B

iter = size(mle_B,4);
mleB = mle_B(:,:,:,iter);
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

theta1 = (2*pi*osc_freq(1))*(1/fs);
mat1 = [cos(theta1) -sin(theta1); sin(theta1) cos(theta1)];
A = blkdiag(rho(1,1)*mat1, rho(2,1)*mat1);

% Loop through each state
for st = 1:numStates
    Bj = mleB(:, :, st);
    Q = var_proc_noi * eye(2*numOsci);
    R = var_obs_noi * eye(numNodes);
    
    % Convert frequency to angular frequency in radians/sample
    omega = 2*pi*f/fs; 

    Hyy_f = zeros(numNodes, numNodes, length(f));
    Hy_f = zeros(numNodes, length(f));
    Hxx_f = zeros(2*numOsci, 2*numOsci, length(f));
    Hx_f = zeros(2*numOsci, length(f));

    for idx = 1:length(omega)
        inv_term = inv(eye(size(A)) - A * exp(-1i*omega(idx)));
        Hxx_f(:, :, idx) = (1/fs) * inv_term * Q * inv_term';
        Hyy_f(:, :, idx) = Bj * Hxx_f(:, :, idx) * Bj' + (1/fs) * R;
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
            S_ij = S_y1y2(i, j, freqIndex, st);
            P_ii = S_y(i, freqIndex, st);
            P_jj = S_y(j, freqIndex, st);
            
            % Compute coherence
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
    title(sprintf('B_%g',i))
end
sgtitle('Coherence (MLE B)')


%% Theoretical coherogram for MLE B

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
coherogram_E12_mleB = zeros(N, length(f));
% weighted coherence
for i = 1:N
    coherogram_E12_mleB(i, :) = coherence_E12*sw(i,:)';
end

% Check if all imaginary parts are essentially zero
imagParts = imag(coherogram_E12_mleB);
areAllImagPartsZero = all(abs(imagParts(:)) < 1e-10); 
if areAllImagPartsZero
    disp('All complex parts are zero.');
else
    disp('There are non-zero complex parts.');
end

% plot the coherogram 
figure
plot_matrix(real(coherogram_E12_mleB), time, f, 'n');
colorbar;
axis xy;
%caxis([0 1]); 
xlabel('Time [s]');
ylabel('Frequency [Hz]');
title('Coherogram (MLE B)')


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

upp = max([max(real(coherogram_E12_true)), max(real(coherogram_E12_mleB))])+0.1;
fig = figure;
set(fig, 'Position', [100, 100, 800, 650]);
ax1 = subplot(4,1,1);
plot_matrix(real(coherogram_E12_true), time, f, 'n');
colorbar;
axis xy;
caxis([0 upp]); 
ylabel('True B', 'FontSize',15);xlabel(''); title('');
text('Units', 'normalized', 'Position', [-0.12, 0.95], 'String', 'B', ...
    'FontSize', 30, 'FontWeight', 'bold');

ax2 = subplot(4,1,2);
plot_matrix(C12, t, f, 'n');
colorbar;
axis xy;
caxis([0 1]); 
ylabel('Multitaper', 'FontSize',15);xlabel(''); title('');

ax3 = subplot(4,1,3);
plot_matrix(real(coherogram_E12_mleB), time, f, 'n');
colorbar;
axis xy;
caxis([0 upp]); 
xlabel('Time [s]'); title('');
ylabel('Estimated B', 'FontSize',15);
text(-0.1, 1.88, 'Frequency [Hz]', 'Units', 'normalized', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
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


%% Plot MLE B - magnitude & phase

clear mag; clear phase;
[mag(:,:,1), phase(:,:,:,1)] = funcs2.B_mag_phase(mle_B(:,:,reorder(1),iter));
[mag(:,:,2), phase(:,:,:,2)] = funcs2.B_mag_phase(mle_B(:,:,reorder(2),iter));
[mag(:,:,3), phase(:,:,:,3)] = funcs2.B_mag_phase(mle_B(:,:,reorder(3),iter));

% make magnitude as a 2x2 grid
out = zeros(2,2,2,3);
for k=1:3
    for i=1:2
        out(:,:,i,k) = reshape(mag(:,i,k),2,2);
    end
end

upp = max(out,[],'all'); 
nm = ['1','2','3','4']; 
figure
tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact'); 
for k = 1:3
    for i = 1:2
        nexttile
        imagesc(out(:,:,i,k)); 
        axis off 
        colorbar 
        caxis([0 0.4]); 
        yticks(1:n); xticks(1:n); 
        set(gca, 'YTickLabel', 1:n); 
        set(gca, 'XTickLabel', 1:n); 
        r = 1; 
        for p = 1:2
            for j = 1:2
                text(p, j, nm(r), 'HorizontalAlignment', 'center');
                r = r + 1;
            end
        end
        if k == 1
            title(sprintf('Oscillator %g', i));
        end
        if i == 1
            h=ylabel(sprintf('B_%g', k),'Rotation', 0, 'FontWeight', 'bold');
            currentPosition = get(h, 'Position');
            newPosition = currentPosition + [-0.8, 0, 0]; 
            set(h, 'Position', newPosition);
        end
    end
end
sgtitle('Estimated B - Magnitude')

% phase
figure
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
for k = 1:3
    for i = 1:2
        nexttile
        funcs2.imagesc2(phase(:,:,i,k))
        axis on
        yticks(1:n); xticks(1:n); 
        set(gca, 'YTickLabel', 1:n); 
        set(gca, 'XTickLabel', 1:n); 
        colorbar
        colormap(hsv)
        caxis([-pi pi])
        if k == 1
            title(sprintf('Oscillator %g', i));
        end
    end
end
sgtitle('Estimated B - Phase')


%% Plot true B - magnitude & phase

clear true_mag; clear true_phase;
[true_mag(:,:,1), true_phase(:,:,:,1)] = funcs2.B_mag_phase(squeeze(B(1,:,:)));
[true_mag(:,:,2), true_phase(:,:,:,2)] = funcs2.B_mag_phase(squeeze(B(2,:,:)));
[true_mag(:,:,3), true_phase(:,:,:,3)] = funcs2.B_mag_phase(squeeze(B(3,:,:)));

true_out = zeros(2,2,2,3);
for k=1:3
    for i=1:2
        true_out(:,:,i,k) = reshape(true_mag(:,i,k),2,2);
    end
end
true_upp = max(true_out,[],'all');

figure
tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact'); 
for k = 1:3
    for i = 1:2
        nexttile
        imagesc(true_out(:,:,i,k)); 
        axis off 
        colorbar 
        caxis([0 true_upp]); 
        yticks(1:n); xticks(1:n); 
        set(gca, 'YTickLabel', 1:n); 
        set(gca, 'XTickLabel', 1:n); 
        r = 1; 
        for p = 1:2
            for j = 1:2
                text(p, j, nm(r), 'HorizontalAlignment', 'center');
                r = r + 1;
            end
        end
        if k == 1
            title(sprintf('Oscillator %g', i));
        end
    end
end
sgtitle('True B - Magnitude')

% phase
figure
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
for k = 1:3
    for i = 1:2
        nexttile
        funcs2.imagesc2(true_phase(:,:,i,k))
        axis on
        yticks(1:n); xticks(1:n); 
        set(gca, 'YTickLabel', 1:n); 
        set(gca, 'XTickLabel', 1:n); 
        colorbar
        colormap(hsv)
        caxis([-pi pi])
        if k == 1
            title(sprintf('Oscillator %g', i));
        end 
    end
end
sgtitle('True B - Phase')


%% Plot mag & phase for true A and MLE A together
upp_true = max(true_mag,[],'all');
upp = max(mag,[],'all');

fig = figure;
set(fig, 'Position', [100, 100, 600, 650]);
tiledlayout(6,4); 
for i=1:numStates
    for j=1:2
        r=1;
        nexttile
        imagesc(true_out(:,:,j,i))
        ax = gca;
        ax.XTick = [];
        ax.YTick = [];
        for q=1:2
            for w=1:2
                text(q,w, sprintf('%g',r),'HorizontalAlignment', 'center')
                r = r+1;
            end
        end
        caxis([0 upp_true])
        colormap(gca, 'parula')
        colorbar
        if i==1
            if j==1
                title('Oscillator 1')
                ylabel('State 1')
                colorbar('off')
            else
                title('Oscillator 2')
            end
        elseif i==2 && j==1
            ylabel('State 2')
            colorbar('off')
        elseif i==3 && j==1
            ylabel('State 3')
            colorbar('off')
        end
    end
    for j=1:2
        r=1;
        nexttile
        imagesc(out(:,:,j,i))
        ax = gca;
        ax.XTick = [];
        ax.YTick = [];
        for q=1:2
            for w=1:2
                text(q,w, sprintf('%g',r),'HorizontalAlignment', 'center')
                r = r+1;
            end
        end
        colorbar
        caxis([0 upp])
        colormap(gca, 'parula')
        if j==1
            colorbar('off')
        end
        if i==1 && j==1
            title('Oscillator 1')
        elseif i==1 && j==2
            title('Oscillator 2')
        end
    end
end
text('Units', 'normalized', 'Position', [-6.2, 4], 'String', 'D', ...
     'FontSize', 30, 'FontWeight', 'bold');
text(-4.1, 4, 'True B', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'FontSize', 15,...
     'VerticalAlignment', 'bottom');
text(-0.3, 4, 'Estimated B', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'FontSize', 15,...
     'VerticalAlignment', 'bottom');
text(-6.15, 1.8, 'Magnitude', 'Units', 'normalized', 'Rotation', 90, ...
     'HorizontalAlignment', 'center', 'FontSize', 15,...
     'VerticalAlignment', 'bottom');
for i=1:numStates
    for j = 1:2
        nexttile
        funcs2.imagesc2(true_phase(:,:,j,i))
        axis on
        yticks(1:n); xticks(1:n);
        set(gca, 'YTickLabel', 1:n);
        set(gca, 'XTickLabel', 1:n);
        colormap(gca, 'hsv')
        caxis([-pi pi])
        cb = colorbar;
        cb.Ticks = [-pi 0 pi];
        cb.TickLabels = {'-π', '0', 'π'};
        if i==1 && j==1
            ylabel('State 1')
            colorbar('off')
        elseif i==2 && j==1
            ylabel('State 2')
            colorbar('off')
        elseif i==3 && j==1
            ylabel('State 3')
            colorbar('off')
        end
    end
    for j = 1:2
        nexttile
        funcs2.imagesc2(phase(:,:,j,i))
        axis on
        yticks(1:n); xticks(1:n);
        set(gca, 'YTickLabel', 1:n);
        set(gca, 'XTickLabel', 1:n);
        colormap(gca, 'hsv')
        caxis([-pi pi])
        cb = colorbar;
        cb.Ticks = [-pi 0 pi];
        cb.TickLabels = {'-π', '0', 'π'};
        if j==1
            colorbar('off')
        end
    end
end
text('Units', 'normalized', 'Position', [-6.6, 3.7], 'String', 'E', ...
     'FontSize', 30, 'FontWeight', 'bold');
text(-6.35, 1.86, 'Phase', 'Units', 'normalized', 'Rotation', 90, ...
     'HorizontalAlignment', 'center', 'FontSize', 15,...
     'VerticalAlignment', 'bottom');


