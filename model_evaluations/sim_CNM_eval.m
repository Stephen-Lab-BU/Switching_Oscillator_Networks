% Correlated Noise Model
% ========================================================================= 
% Author:          Wan-Chi Hsin
% Date:            April 16, 2024
% MATLAB Version:  R2021b
% =========================================================================
% simulation study
% single rhythm
% 10 electrodes
% 3 switching states
% random network structure in each state
% network structure is reflected through the process noise cov matrix (Q/Sigma)

%% Setting
rng(35)
addpath(genpath(pwd))

% ===== tuning parameter =====
var_obs_noi = 8;        %var(observation noise)
aAR = 0.8;              %AR asso. w/ A
var_proc_noi = 1;       %var(process noise) on the diagonal

numNodes = 10;          %number of electrodes
numOsci = numNodes;     %number of oscillators
numStates = 3;          %number of switching states

fs = 100;               %sampling frequency
osciFreq = 7;           %oscillation frequency
T = 300; TT=T;          %total time
zij = 0.0001/2;         %prob of switching 
percEdges = [0.05 0.05 0.05];   % Perc of edges for each state
mag_range = [0.6 0.4];          % range of linkage strength


% ===== variables setup =====
dt = 1/fs;              %sampling interval
N = T*fs;               %number of points in data

% ** process transition matrix: A matrix **
A = zeros(2*numOsci,2*numOsci);
w = 2*pi*osciFreq/fs;
for i = 1:numOsci
    block = aAR*[cos(w), -sin(w); sin(w), cos(w)];
    A(2*i-1:2*i, 2*i-1:2*i) = block;
end

% ** observation matrix: B matrix **
B = zeros(numNodes,2*numOsci);
for i=1:numOsci
    B(i,2*i-1) = 1;
end

% ** transition matrix Z **
Z = [1-2*zij zij zij; zij 1-2*zij zij; zij zij 1-2*zij];

% ** switching variable S **
S_t = zeros(1, N);
S_t(1) = randi([1, 3], 1, 1);   % Initialize the starting state
for t = 2:N
    currentState = S_t(t-1);
    nextState = randsample(1:3, 1, true, Z(currentState, :));
    S_t(t) = nextState;
end
S_t = S_t';

% ===== process noise cov matrix: Q/Sigma matrix =====
rng(22)
% Define adjacency matrix for linkages
adjMatrix = zeros(numNodes, numNodes, numStates);

% Calculate the total number of edges
totalEdges = round(percEdges * numNodes * (numNodes - 1));

for i=1:numStates
    % Randomly assign edges 
    edgesAdded = 0;
    while edgesAdded < totalEdges(i)
        % Randomly select a source and target node
        source = randi([1, numNodes]);
        target = randi([1, numNodes]);
        
        % Ensure source is not equal to target and no edge already exists
        if source ~= target && adjMatrix(source, target,i) == 0 && adjMatrix(target, source,i) ~= 1
            adjMatrix(source, target,i) = 1; % Add edge
            edgesAdded = edgesAdded + 1;
        end
    end 
end
adjMatrix = updateAdjMatrix(adjMatrix, numNodes, numStates);

% Define magnitude and phase for the linkages
magMatrix = zeros(size(adjMatrix));
phaseMatrix = zeros(size(adjMatrix));
mag_upper = mag_range(1);
mag_lower = mag_range(2);
commonPhase = 0;

for j = 1:numStates
    % Find indices where adjMatrix is 1 to determine the number of edges in the state
    [edgeRows, edgeCols] = find(adjMatrix(:,:,j) == 1);
    numEdges = length(edgeRows); % Number of edges in this state
    
    % Assign mag and phase values to edges symmetrically
    for i = 1:numEdges
        % Generate a magnitude value for the current edge
        currentMag = mag_upper - (mag_upper - mag_lower) * rand();

        % Assign the magnitude value symmetrically for both directions of the edge
        magMatrix(edgeRows(i), edgeCols(i), j) = currentMag;
        magMatrix(edgeCols(i), edgeRows(i), j) = currentMag; 

        % Assign the phase value
        phaseMatrix(edgeRows(i), edgeCols(i), j) = commonPhase;
    end
end


% Define Q matrix
QMatrix = zeros(numNodes*2, numNodes*2,numStates);
% Iterate through each element in the magMatrix and phaseMatrix
for k=1:numStates
    for i = 1:size(magMatrix, 1)
        for j = 1:size(magMatrix, 2)
            % Extract mag_val and phase_val for the current edge
            mag_val = magMatrix(i, j,k);
            phase_val = phaseMatrix(i, j,k);
            
            % Calculate the 2x2 block using mag_val and phase_val
            block = mag_val * [cos(phase_val) -sin(phase_val); sin(phase_val) cos(phase_val)];
            
            % Calculate the position in the combinedMatrix
            rowStart = (i-1)*2 + 1;
            colStart = (j-1)*2 + 1;
            
            % Place the block in the correct position in the combinedMatrix
            QMatrix(rowStart:rowStart+1, colStart:colStart+1,k) = block;
        end
        QMatrix(rowStart:rowStart+1, rowStart:rowStart+1,k) = var_proc_noi*eye(2);
    end
end


% ===== Simulate the system =====
x = zeros(2*numOsci, N);
y = zeros(numNodes,N);

% Initialize oscillators on the unit circle using the random angles
initial_phase =  2*pi*rand(1, numOsci);
cos_vals = cos(initial_phase);
sin_vals = sin(initial_phase);
for i=1:numOsci
    x(2*i-1,1) = cos_vals(i);
    x(2*i,1) = sin_vals(i);
end

% Simulate the system
v_t = zeros(numNodes,N);
for t = 2:N
    currentState = S_t(t); % Get the current state from S_t to decide which Q matrix to use
    QMatrixCurrent = QMatrix(:,:,currentState);
    % Ensure QMatrixCurrent is symmetric and positive definite for Cholesky decomposition
    % Use 'lower' since chol returns lower triangular matrix by default in MATLAB
    QChol = chol(QMatrixCurrent, 'lower');
    
    % Generate process noise u_t using the Cholesky decomposition of the current Q matrix
    u_t = QChol * randn(2*numOsci, 1); 

    v_t(:,t) = sqrt(var_obs_noi) * randn(numNodes, 1); % Observation noise
    x(:, t) = A * x(:, t-1) + u_t; % Use A matrix corresponding to current state
    y(:, t) = B * x(:, t) + v_t(:,t); % Observation equation
end


%% Coherogram - multitaper
% *** make sure you have chronux installed ***
% only show the pair of 1 & 2 here
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
[C,Cphi,Ccs,C1,C2,t,f] = cohgramc(data(1,:)', data(2,:)', movingwin, params); 


%% Theoretical cross-spectrum & coherence for true Q

% Initialize the coherence matrix for all states
Coh = zeros(numNodes, numNodes, numStates);
S_y = zeros(numNodes, length(f), numStates);
S_y1y2 = zeros(numNodes, numNodes,length(f), numStates);

% Loop through each state
for st = 1:numStates
    Qj = QMatrix(:, :, st);
    R = var_obs_noi * eye(numNodes);

    omega = 2*pi*f/fs; % Convert frequency to angular frequency in radians/sample

    Hyy_f = zeros(numNodes, numNodes, length(f));
    Hy_f = zeros(numNodes, length(f));
    Hxx_f = zeros(2*numOsci, 2*numOsci, length(f));
    Hx_f = zeros(2*numOsci, length(f));

    for idx = 1:length(omega)
        inv_term = inv(eye(size(A)) - A * exp(-1i*omega(idx)));
        Hxx_f(:, :, idx) = (1/fs) * inv_term * Qj * inv_term';
        Hyy_f(:, :, idx) = B * Hxx_f(:, :, idx) * B' + (1/fs) * R;
        Hx_f(:, idx) = diag(Hxx_f(:, :, idx));
        Hy_f(:, idx) = diag(Hyy_f(:, :, idx));
    end

    % Compute theoretical coherence
    [~, freqIndex] = min(abs(f - osciFreq));

    S_y1y2(:,:,:,st) = Hyy_f;
    S_y(:,:,st) = Hy_f; % Power spectra for each electrode (example data)
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
true_Coh = Coh;
S_y1y2_7hz = squeeze(S_y1y2(:,:,freqIndex,:));


%plot coherence matrix - True Q
figure
tiledlayout(1,numStates)
for i=1:numStates
    nexttile
    n = numNodes;
    imagesc(Coh(:,:,i))
    yticks(1:n); xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colorbar
    caxis([0 1])
    colormap(parula)
    title(sprintf('\\Sigma%g',i))
end
sgtitle(sprintf('Coherence: True \\Sigma (\\alpha=%g, var(R)=%g)', aAR, var_obs_noi))

%plot cross-spectrum matrix - True Q
figure
tiledlayout(1,numStates)
for i=1:numStates
    nexttile
    n = numNodes;
    imagesc(abs(S_y1y2_7hz(:,:,i)))
    yticks(1:n); xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colorbar
    %caxis([0 1])
    colormap(parula)
    title(sprintf('\\Sigma%g',i))
end
sgtitle(sprintf('Cross-spectrum: True \\Sigma (\\alpha=%g, var(R)=%g)', aAR, var_obs_noi))


%% Windowed multitaper coherence -- test if link is significant
alpha = 0.05;

W = 1;
T = 20;
TW = T*W;                        %time-bandwidth product given 2W=2, 
ntapers = 2*TW-1;                %... 2W is the desired freq resolution
params.Fs = 1/dt;                %... 19 tapers are used
params.tapers = [TW,ntapers];    
params.pad = -1;       
params.err=[2 0.05];
window_size = 1;
window_step = 1;
movingwin=[window_size window_step];   %set the moving window dimensions: [window winstep]

[C,Cphi,S12,S1,S2,t,f] = cohgramc(y(1,:)', y(2,:)', movingwin, params); 

numWindows = TT/window_size;
C = zeros(numNodes, numNodes, numWindows, length(f)); 
S12 = zeros(numNodes, numNodes, numWindows, length(f)); 
S1 = zeros(numNodes, numWindows, length(f)); 
Cxy_at_7Hz = zeros(numNodes, numNodes, numWindows);
S12_at_7Hz = zeros(numNodes, numNodes, numWindows);
S1_at_7Hz = zeros(numNodes, numWindows); 
F_statistic = zeros(numNodes, numNodes, numWindows);
p_value = zeros(numNodes, numNodes, numWindows);
test_output = zeros(numNodes, numNodes, numWindows);

n = numNodes;
for i = 1:n % Loop over all electrodes except the last one
    for j = 1:n % Loop over electrodes after the current i to ensure j > i
        [C(i,j,:,:),Cphi,S12(i,j,:,:),S1(i,:,:),S2,t,f] = cohgramc(y(i,:)', y(j,:)', movingwin, params); 

        % Extract the coherence value at 7 Hz
        [~, idx] = min(abs(f - osciFreq));
        Cxy_at_7Hz = C(:,:,:,idx);
        S12_at_7Hz = S12(:,:,:,idx);
        S1_at_7Hz = S1(:,:,idx);

        for k=1:numWindows
            % Extract the coherence value at 7 Hz for each window
            windowed_val = Cxy_at_7Hz(i,j,k);
            % Compute the F-test statistic
            F_statistic(i,j,k) = ntapers * windowed_val^2/(1-windowed_val^2);
            % Compute the p-value
            p_value(i,j,k) = 1 - fcdf(F_statistic(i,j,k),2,2*ntapers);

            alpha_adjusted = alpha;
            if p_value(i,j,k) < alpha_adjusted
                test_output(i,j,k) = 1;
            end
        end
    end
end

% ======== testing coherence for each window ========
% F-test
% create true network for each window 
S_tWindowed = zeros(numWindows,1);
window_samples = window_size * fs;

% Process each window
true_output = zeros(numNodes,numNodes,numWindows);
for i = 1:numWindows
    start_index = (i-1) * window_samples + 1;
    end_index = start_index + window_samples - 1;
    
    % Handle case where S_t does not perfectly divide into windows
    if end_index > length(S_t)
        end_index = length(S_t);
    end
    % Extract the window
    window = S_t(start_index:end_index);
    
    % Determine the most frequent value in the window
    values = [1, 2, 3];
    counts = arrayfun(@(x) sum(window == x), values);
    [~, idx] = max(counts);
    most_frequent_value = values(idx);
    
    % Assign the most frequent value to the S_tWindowed
    S_tWindowed(i) = most_frequent_value;
    true_output(:,:,i) = true_Coh(:,:,S_tWindowed(i));
    %true_output(:,:,i) = triu(true_Coh(:,:,S_tWindowed(i))) - eye(numNodes);
end
true_output(true_output > 0.1) = 1;

TP = 0; FP = 0;
for i=1:numWindows
    % True Positives
    TP = TP + sum((true_output(:,:,i) == 1) & (test_output(:,:,i) == 1),'all') - numNodes; % do not count diagonal
    % False Positives
    FP = FP + sum((true_output(:,:,i) == 0) & (test_output(:,:,i) == 1),'all'); 
end
totalActualPos = sum(true_output == 1,'all') - numNodes*numWindows; 
totalActualNeg = sum(true_output ~= 1,'all'); 
fracCorrect = TP / totalActualPos; % Fraction of correctly identified 1's
falsePosRate = FP / totalActualNeg; % False Positive Rate

fprintf('Sensitivity: %.4f (%g/%g)\n', fracCorrect, TP, totalActualPos);
fprintf('False Positive Rate: %.4f (%g/%g)\n', falsePosRate, FP,totalActualNeg);


%% Norm (multitaper vs true cross-spectrum) 
diff = zeros(numNodes,numNodes,numWindows);
diag_norm = zeros(numWindows,1);
offdiag_norm = zeros(numWindows,1);
mask = ~eye(numNodes, numNodes);

% S12_at_7Hz: multitaper 
% S_y1y2_7hz: true cross-spectrum, which one to use depending on S_tWindowed

for i=1:numWindows
    % for each window, compute the diff btw the true CS and the MT SC
    diff(:,:,i) =  S12_at_7Hz(:,:,i) - S_y1y2_7hz(:,:,S_tWindowed(i));
    slice = diff(:,:,i);
    % diagonal: norm
    diag_val = diag(slice);
    diag_norm(i) = norm(abs(diag_val))/sqrt(length(diag_val));
    % off-diagonal: Frobenius norm
    offdiag_val = slice(mask);
    offdiag_norm(i) = norm(abs(offdiag_val), 'fro')/sqrt(numel(offdiag_val));
end
meanDiag = mean(diag_norm);
sdDiag = std(diag_norm);
meanOffdiag = mean(offdiag_norm);
sdOffdiag = std(offdiag_norm);

fprintf('The Frobenius norm is %0.4f (%0.4f) (diagonal), %0.4f (%0.4f) (off-diagonal).\n', ...
        meanDiag, sdDiag, meanOffdiag, sdOffdiag);


%% EM on Q
rng(2238)

iter = 30; 
n = numNodes;
k = numOsci;
M = numStates;
Q = zeros(2*k,2*k,M);
R = zeros(n,n,M);

% Initial Q matrix 
Q(:,:,1) = funcs2.randq(k);          
Q(:,:,2) = funcs2.randq(k);
Q(:,:,3) = funcs2.randq(k);

% A matrix
clear A0;
A0(:,:,1) = A;   
A0(:,:,2) = A; 
A0(:,:,3) = A;

% B matrix
B0 = zeros(n,2*k,M);
B0(:,:,1) = B;   
B0(:,:,2) = B; 
B0(:,:,3) = B;

% R matrix
R(:,:,1) = var_obs_noi*eye(n);
R(:,:,2) = var_obs_noi*eye(n);
R(:,:,3) = var_obs_noi*eye(n);

X_0 = mvnrnd(zeros(2*k,1), eye(2*k))';
C = Z;  %transition matrix

% EM on Q
[mle_Q,X_RTS,SW,Q_func] = em_projQ(y',1e-5,iter,A0,B0,Q,R,C,X_0);


%% EM on A
rng(2238)

numOsci = numNodes;

iter = 20; 
n = numNodes;
k = numOsci;
M = numStates;
A = zeros(2*k,2*k,M);
B = zeros(n,2*k,M);
Q = zeros(2*k,2*k,M);
R = zeros(n,n,M);

% Initial A matrix
A(:,:,1) = funcs2.randa2(k,osciFreq,fs,aAR);      
A(:,:,2) = funcs2.randa2(k,osciFreq,fs,aAR);    
A(:,:,3) = funcs2.randa2(k,osciFreq,fs,aAR);   

% B matrix
for i=1:n
    B(i,2*i-1:2*i,:) = 1/sqrt(2);
end

% R matrix
R(:,:,1) = var_obs_noi*eye(n);
R(:,:,2) = var_obs_noi*eye(n);
R(:,:,3) = var_obs_noi*eye(n);

% Q matrix
Q(:,:,1) = var_proc_noi*eye(2*k);        
Q(:,:,2) = var_proc_noi*eye(2*k);       
Q(:,:,3) = var_proc_noi*eye(2*k);       

X_0 = mvnrnd(zeros(2*k,1), eye(2*k))';
C = Z;  %transition matrix

% EM on A
[mle_A,X_RTS,V_RTS,SW,Q_func] = em_projA(y',1e-5,iter,A,B,Q,R,C,X_0);


%% EM on B
rng(22)

var_obs_noi = 20;
numOsci = 2;
iter = 20;
n = numNodes;
k = numOsci;
M = numStates;
Q = zeros(2*k,2*k,M);
R = zeros(n,n,M);

% Initial B matrix
% everything is randomly weakly correlated
B0 = .05*rand(n,2*k,M);    

% A matrix
clear A0;
A = zeros(2*numOsci,2*numOsci);
w = 2*pi*osciFreq/fs;
for i = 1:numOsci
    block = aAR*[cos(w), -sin(w); sin(w), cos(w)];
    A(2*i-1:2*i, 2*i-1:2*i) = block;
end
A0(:,:,1) = A;
A0(:,:,2) = A;
A0(:,:,3) = A;

% R matrix
R(:,:,1) = var_obs_noi*eye(n);
R(:,:,2) = var_obs_noi*eye(n);
R(:,:,3) = var_obs_noi*eye(n);

% Q matrix
Q(:,:,1) = var_proc_noi*eye(2*k);        
Q(:,:,2) = var_proc_noi*eye(2*k);       
Q(:,:,3) = var_proc_noi*eye(2*k);       
  

X_0 = mvnrnd(zeros(2*k,1), eye(2*k))';
C = Z;  %transition matrix

% EM on B
[mle_B,X_RTS,SW,Q_func] = em_B(y',1e-5,iter,A0,B0,Q,R,C,X_0);



%defined functions
function undirAdjMatrix = updateAdjMatrix(adjMatrix, numNodes, numStates)
    undirAdjMatrix = zeros(numNodes, numNodes, numStates);
    for i = 1:numStates
        % Create a graph representation for each state
        G = graph(adjMatrix(:,:,i) + adjMatrix(:,:,i)');
        for j = 1:numNodes
            for k = (j+1):numNodes % to avoid repetition and self-loops
                if shortestpath(G, j, k, 'Method', 'unweighted') ~= Inf
                    % A path exists between j and k, update undirAdjMatrix accordingly
                    undirAdjMatrix(j, k, i) = 1;
                    undirAdjMatrix(k, j, i) = 1;
                end
            end
        end
    end
end
