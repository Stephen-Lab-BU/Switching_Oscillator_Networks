% Directed Influence Model
% ========================================================================= 
% Author:          Wan-Chi Hsin
% Date:            April 16, 2024
% MATLAB Version:  R2021b
% =========================================================================
% simulation study
% network structure is reflected through the state transition matrix (A)
% use this file for 
% (1) switching accuracy
% (2) the element-wise error between the estimated and true
% cross-spectal matrices (off-diagonal elements)
% (3) tests for coherence between pair of nodes 

% ** first, load results from EM on A **

%% Switching accuracy
t = (dt:dt:TT);
iter = size(mle_A,4);
sw = SW(:,:,iter+1);
threshold = 0.05;

figure
tiledlayout(numStates,1);
for i=1:numStates
    nexttile
    plot(t,sw(:,i),'LineWidth',2)
    set(gca,'YTick',[0 1],'YLim',[-0.1 1.1]);
    title(sprintf('A%g', i));
end
xlabel('Time [s]')
ylabel('Probability')
sgtitle('Estimated moment by moment probability of switching states')

% Switching accuracy
estSt = zeros(N, 1);
estStConfidence = zeros(N, 1); % To track confidence of estimation

% Step 1: Determine the most likely state with confidence threshold
for t=1:N
    [~, estSt(t)] = max(sw(t,:));
end
for t = 1:N
    [probValues, probIndices] = sort(sw(t, :), 'descend');
    if probValues(1) - probValues(2) > threshold
        estSt(t) = probIndices(1);
        estStConfidence(t) = 1; % Confident estimation
    else
        %estSt(t) = 0; % Uncertain estimation
        estStConfidence(t) = 0;
    end
end

% Step 2: Reorder based on the most frequent values where S_t equals i
mapping = zeros(1,3);
for i = 1:3
    modeVal = mode(estSt(S_t == i));
    mapping(i) = modeVal;
end
disp(mapping);
reorder = mapping;

% Step 3: Compute the switching estimates accuracy considering only confident estimates
estStreorder = zeros(size(estSt));
for i = 1:length(estSt)
    if estStConfidence(i) % Consider only confident estimates
        idx = find(reorder == estSt(i));
        if length(idx) == 1
            estStreorder(i) = idx;
        else
            filteredValues = S_t(ismember(S_t, idx));
            estStreorder(i) = mode(filteredValues);
        end
    else
        estStreorder(i) = 0; % Uncertain estimate remains marked as such
    end
end

% Calculate accuracy only for confident estimations
sameValuesLogical = (estStreorder == S_t) & estStConfidence;
fracSWcorrect = mean(sameValuesLogical);
numera = sum(sameValuesLogical);
fprintf('Fraction of accurate, confident estimates: %.3f (%g/%g) \n', fracSWcorrect, numera, N);


%% Theoretical coherence & cross-spectrum of MLE A

% Initialize the coherence matrix for all states
iter = size(mle_A,4);
CohMLEA = zeros(numNodes, numNodes, numStates);
S_y = zeros(numNodes, length(f), numStates);
S_y1y2 = zeros(numNodes, numNodes,length(f), numStates);

numOsci = numNodes;
% B matrix
B = zeros(numNodes,2*numOsci);
for i=1:n
    B(i,2*i-1:2*i) = 1/sqrt(2);
end

% Loop through each state
for st = 1:numStates
    Aj = mle_A(:, :, st, iter);
    Q = var_proc_noi * eye(2*numOsci);
    R = var_obs_noi * eye(numNodes);

    omega = 2*pi*f/fs; % Convert frequency to angular frequency in radians/sample

    Hyy_f = zeros(numNodes, numNodes, length(f));
    Hy_f = zeros(numNodes, length(f));
    Hxx_f = zeros(2*numOsci, 2*numOsci, length(f));
    Hx_f = zeros(2*numOsci, length(f));

    for idx = 1:length(omega)
        inv_term = inv(eye(size(Aj)) - Aj * exp(-1i*omega(idx)));
        Hxx_f(:, :, idx) = (1/fs) * inv_term * Q * inv_term';
        Hyy_f(:, :, idx) = B * Hxx_f(:, :, idx) * B' + (1/fs) * R;
        Hx_f(:, idx) = diag(Hxx_f(:, :, idx));
        Hy_f(:, idx) = diag(Hyy_f(:, :, idx));
    end

    % Compute theoretical coherence
    [~, freqIndex] = min(abs(f - osciFreq));

    S_y1y2(:,:,:,st) = Hyy_f;
    S_y(:,:,st) = Hy_f;% Power spectra for each electrode (example data)
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
    CohMLEA(:, :, st) = abs(coherenceMatrix);
end
S_y1y2_7hz_MLEA = squeeze(S_y1y2(:,:,freqIndex,:));


% plot coherence matrix - MLE A
figure
tiledlayout(1,numStates)
for i=1:numStates
    nexttile
    n = numNodes;
    imagesc(CohMLEA(:,:,reorder(i)))
    yticks(1:n); xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colorbar
    caxis([0 1])
    colormap(parula)
    title(sprintf('A%g',i))
end
sgtitle(sprintf('Coherence: MLE A (\\alpha=%g, var(R)=%g)', aAR, var_obs_noi))


%plot cross-spectrum matrix - MLE A
figure
tiledlayout(1,numStates)
for i=1:numStates
    nexttile
    n = numNodes;
    imagesc(abs(S_y1y2_7hz_MLEA(:,:,reorder(i))))
    yticks(1:n); xticks(1:n);
    set(gca, 'YTickLabel', 1:n)
    set(gca, 'XTickLabel', 1:n)
    colorbar
    %caxis([0 1])
    colormap(parula)
    title(sprintf('A%g',i))
end
sgtitle(sprintf('Cross-spectrum: MLE A (\\alpha=%g, var(R)=%g)', aAR, var_obs_noi))


%% Norm (cross-spectrum)
diff = zeros(numNodes,numNodes,N);
diag_norm = zeros(N,1);
offdiag_norm = zeros(N,1);
mask = ~eye(numNodes, numNodes);

for i=1:N
    wt = sw(i,:);
    weighted_mat = zeros(numNodes,numNodes);
    for j = 1:3
        weighted_mat = weighted_mat + S_y1y2_7hz_MLEA(:,:,j) * wt(j);
    end
    true_mat = S_y1y2_7hz(:,:,S_t(i));
    diff_mat = weighted_mat - true_mat;
    % diagonal - norm
    diag_val = diag(diff_mat);
    diag_norm(i) = norm(abs(diag_val))/sqrt(length(diag_val));
    % off-diagonal - Frobenius norm
    offdiag_val = diff_mat(mask);
    offdiag_norm(i) = norm(abs(offdiag_val), 'fro')/sqrt(numel(offdiag_val));
end
meanDiag = mean(diag_norm);
sdDiag = std(diag_norm);
meanOffdiag = mean(offdiag_norm);
sdOffdiag = std(offdiag_norm);

fprintf('The Frobenius norm (cross-spectrum) is %0.4f (%0.4f) (diagonal), %0.4f (%0.4f) (off-diagonal).\n', ...
        meanDiag, sdDiag, meanOffdiag, sdOffdiag);


%% Test links (coherence matrix)
% fit an empirical gamma distribution

alpha = 0.05;

n = numNodes;
k = numOsci;
M = numStates;
N = TT*fs;

CohMLEArd = zeros(n,n,M);
for i=1:M
    CohMLEArd(:,:,i) = CohMLEA(:,:,reorder(i));
end
maskTrueLinks = true_Coh > 0.05;
maskTrueLinks = double(maskTrueLinks) - eye(n);
maskTrueLinks = logical(maskTrueLinks);
GivenTrue = CohMLEArd(maskTrueLinks);

maskNoLinks = true_Coh < 0.05;
maskNoLinks = double(maskNoLinks);
maskNoLinks = logical(maskNoLinks);
GivenNon = CohMLEArd(maskNoLinks);

% % visualization
% figure;
% hold on; 
% binsGivenTrue = 40;
% histogram(GivenTrue, binsGivenTrue, 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.5);
% binsGivenNon = 50; 
% histogram(GivenNon, binsGivenNon, 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.5);
% hold off;
% legend('GivenTrue', 'GivenNon');  
% xlabel('Value'); 
% ylabel('Probability Density');  
% title('Comparison of Distributions (Coherence)'); 

% figure
% subplot(2,1,1)
% histogram(GivenTrue, binsGivenTrue, 'Normalization', 'pdf', 'FaceColor', 'blue', 'FaceAlpha', 0.5);
% subplot(2,1,2)
% histogram(GivenNon, binsGivenNon, 'Normalization', 'pdf', 'FaceColor', 'red', 'FaceAlpha', 0.5);


% fit gamma
maskAll = ones(n,n,M);
for i = 1:M
    maskAll(:,:,i) = maskAll(:,:,i) - eye(n);
end
maskAll = logical(maskAll);
allEdges = CohMLEArd(maskAll);

pd = fitdist(allEdges, 'Gamma');
disp(['Shape (k) parameter: ', num2str(pd.a)]);
disp(['Scale (θ) parameter: ', num2str(pd.b)]);

% Find the critical value for an upper-tailed test
critical_value_upper = gaminv(1 - alpha, pd.a, pd.b);
disp(['Critical Value (upper-tailed): ', num2str(critical_value_upper)]);

% ======== Sensitivity and FPR ========
statistics = CohMLEArd;
test_output2 = zeros(n,n,M);
for m=1:M
    for i=1:n
        for j=1:n
            if i~=j
                if statistics(i,j,m) > critical_value_upper
                    test_output2(i,j,m) = 1;
                end
            end
        end
    end
end

true_links = zeros(numNodes,numNodes,numStates);
for i=1:numStates
    true_links(:,:,i) = real(true_Coh(:,:,i)) - eye(numNodes);
end
TP = sum((test_output2 == 1) & (true_links > 0.01),'all'); % True Positives
FP = sum((test_output2 == 1) & (true_links < 0.01),'all'); % False Positives

totalActualPos = sum(true_links > 0.01,'all'); % Total actual positives in true_val
totalActualNeg = sum(true_links < 0.01,'all')-numNodes*numStates; % Total actual negatives

fracCorrect = TP / totalActualPos; % Fraction of correctly identified 1's
falsePosRate = FP / totalActualNeg; % False Positive Rate

fprintf('Sensitivity: %.3f (%g/%g)\n', fracCorrect, TP, totalActualPos);
fprintf('False Positive Rate: %.3f (%g/%g)\n', falsePosRate, FP,totalActualNeg);

% %% histogram of all edges
% figure
% bins = 100;
% 
% histogram(allEdges, bins, 'Normalization', 'pdf');
% hold on;  
% x_values = linspace(min(allEdges), max(allEdges), 1000);  
% pdf_values = pdf(pd, x_values);
% 
% plot(x_values, pdf_values, 'r-', 'LineWidth', 2); 
% xlabel('Data Values');
% ylabel('Probability Density');
% title('Histogram and Fitted Gamma Distribution');
% legend('Data Histogram', 'Fitted Gamma PDF');
% hold off; 


%% Test links (coherence matrix) - truncated

upperBound = 0.09; 
truncatedData = allEdges(allEdges < upperBound);
pd_truncated = fitdist(truncatedData, 'Gamma');

disp(['Shape parameter (k): ', num2str(pd_truncated.a)]);
disp(['Scale parameter (θ): ', num2str(pd_truncated.b)]);

% Find the critical value for an upper-tailed test
critical_value_upper = gaminv(1 - alpha, pd_truncated.a, pd_truncated.b);
disp(['Critical Value (upper-tailed): ', num2str(critical_value_upper)]);

% ======== Sensitivity and FPR ========
statistics = CohMLEArd;
test_output2 = zeros(n,n,M);
for m=1:M
    for i=1:n
        for j=1:n
            if i~=j
                if statistics(i,j,m) > critical_value_upper
                    test_output2(i,j,m) = 1;
                end
            end
        end
    end
end

true_links = zeros(numNodes,numNodes,numStates);
for i=1:numStates
    true_links(:,:,i) = real(true_Coh(:,:,i)) - eye(numNodes);
end
TP = sum((test_output2 == 1) & (true_links > 0.01),'all'); % True Positives
FP = sum((test_output2 == 1) & (true_links < 0.01),'all'); % False Positives

totalActualPos = sum(true_links > 0.01,'all'); % Total actual positives in true_val
totalActualNeg = sum(true_links < 0.01,'all')-numNodes*numStates; % Total actual negatives

fracCorrect = TP / totalActualPos; % Fraction of correctly identified 1's
falsePosRate = FP / totalActualNeg; % False Positive Rate

fprintf('Sensitivity: %.3f (%g/%g)\n', fracCorrect, TP, totalActualPos);
fprintf('False Positive Rate: %.3f (%g/%g)\n', falsePosRate, FP,totalActualNeg);


varName = ['sen_dim_var', num2str(var_obs_noi)];
eval([varName ' = fracCorrect;']);

fileName = [varName '.mat'];
eval(['save(fileName, ''' varName ''');']);

% %plot histograms
% figure
% bins = 40;
% histogram(allEdges, bins, 'Normalization', 'pdf');
% hold on;  
% x_values = linspace(min(allEdges), max(allEdges), 1000);  
% pdf_values = pdf(pd_truncated, x_values);
% plot(x_values, pdf_values, 'r-', 'LineWidth', 2); 
% xlabel('Data Values');
% ylabel('Probability Density');
% title('Histogram and Fitted Gamma Distribution (Coherence)');
% legend('Data Histogram', 'Fitted Gamma PDF');
% hold off; 


