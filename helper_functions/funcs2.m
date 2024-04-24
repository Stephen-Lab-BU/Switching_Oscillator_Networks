classdef funcs2
    methods(Static)
       
        % ===== create A matrices, 3 states =====
        function A_sim = create_A3(freq, fs, rho, M)
            k = size(freq,2);
            A_sim = zeros(k,M,2,2); % M: # of state
            
            for i=1:size(freq,2)
                theta = (2*pi*freq(i))*(1/fs);
                rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                for j=1:M
                    A_sim(i,j,:,:) = rho(i,j) * rot;
                end
            end
        end

        % ===== create A matrices, 3 states for DIM =====
        function A_sim = create_A3dim(freq,Fs,M,alpha,theta,beta,row_sum,aAR)
            k = size(freq,2);
            A_sim = zeros(2*k,2*k,M); % M: # of state
            for m=1:M
                for i=1:k
                    for j=1:k
                        if i~=j
                            Rmat = [cos(theta(i,j,m)) -sin(theta(i,j,m)); sin(theta(i,j,m)) cos(theta(i,j,m))];
                            A_sim(2*i-1:2*i,2*j-1:2*j,m) = alpha(i,j,m)*beta(i,j,m)*Rmat;
                        else
                            Rmat = [cos(2*pi*freq(i)/Fs) -sin(2*pi*freq(i)/Fs);sin(2*pi*freq(i)/Fs) cos(2*pi*freq(i)/Fs)];
                            A_sim(2*i-1:2*i,2*j-1:2*j,m) = aAR*Rmat - row_sum(i,m)*eye(2);
                        end
                    end
                end
            end
        end
        
        % ===== create carrier signals, x, 3 states for DIM =====
        function x = carrier_signal3dim(k,N,initial_phase,jumppoint1,jumppoint2,sub_A,noi)
            x = zeros(2*k,N);
            x(:,1) = repmat([1,0],[1,k]);
            for t=1:N
                if t <= jumppoint1
                    x(:,t+1) = sub_A(:,:,1) * x(:,t) + noi(:,t,1);
                elseif (t > jumppoint1 && t <= jumppoint2)
                    x(:,t+1) = sub_A(:,:,2) * x(:,t) + noi(:,t,2);
                else
                    x(:,t+1) = sub_A(:,:,3) * x(:,t) + noi(:,t,3);
                end
            end
            x = x(:,2:N+1);
        end

        % ===== create carrier signals, x, 3 states =====
        function x = carrier_signal3(k,N,initial_phase,jumppoint1,jumppoint2,sub_A,noi)
            x = zeros(2*k,N);
            for i=1:k
                % set starting point
                x(i*2-1:i*2,1) = [1, 0];
                for t=1:N
                    if t <= jumppoint1
                        x(i*2-1:i*2,t+1) = squeeze(sub_A(i,1,:,:)) * x(i*2-1:i*2,t) + noi(i*2-1:i*2,t,1);
                    elseif (t > jumppoint1 && t <= jumppoint2)
                        x(i*2-1:i*2,t+1) = squeeze(sub_A(i,2,:,:)) * x(i*2-1:i*2,t) + noi(i*2-1:i*2,t,2);
                    else
                        x(i*2-1:i*2,t+1) = squeeze(sub_A(i,3,:,:)) * x(i*2-1:i*2,t) + noi(i*2-1:i*2,t,3);
                    end 
                end
            end
            
            x = x(:,2:N+1);
        end
        
        % ===== create signals, y, 3 states =====
        function y = obs_signal3(x,B,n,N,jumppoint1,jumppoint2, obs_noi)
            y = zeros(n,N);
            for t=1:N
                if t <= jumppoint1
                    y(:,t) = squeeze(B(1,:,:)) * x(:,t);
                elseif (t > jumppoint1 && t <= jumppoint2)
                    y(:,t) = squeeze(B(2,:,:)) * x(:,t);
                else
                    y(:,t) = squeeze(B(3,:,:)) * x(:,t);
                end
            end
            y = y + obs_noi;
        end
                
        % ===== create random Q matrix =====
        % randomly weakly correlated
        % random phase between -pi to pi
        % random rho between interval (a,b) using a + (b-a)*rand(1)
        function q = randq(k)
            q = eye(2*k);
            a = 0.01;
            b = 0.2;
            for i=1:k
                for j=1:k
                    if i < j
                        rdrho = a + (b-a)*rand(1);
                        rdtheta = -pi + (2*pi)*rand(1);
                      
                        q(i*2-1,j*2-1) = rdrho * cos(rdtheta);
                        q(i*2,j*2-1) = rdrho * sin(rdtheta);
                        q(i*2-1,j*2) = rdrho * -sin(rdtheta);
                        q(i*2,j*2) = rdrho * cos(rdtheta);
                    end
                end
            end
            
            for i=1:k
                for j=1:k
                    if i > j
                        q(i*2-1,j*2-1) = q(j*2-1,i*2-1);
                        q(i*2,j*2-1) = -q(j*2,i*2-1);
                        q(i*2-1,j*2) = -q(j*2-1,i*2);
                        q(i*2,j*2) = q(j*2,i*2);
                    end
                end
            end
        end

        % ===== create random A matrix (aAR=0.9) =====
        % randomly weakly correlated
        % random phase between -pi to pi
        % random rho between interval (a,b) using a + (b-a)*rand(1)
        function a0 = randa(k,osc_freq,fs)
            a0 = eye(2*k);
            a = 0.01;
            b = 0.2;
            % on diagonal - rotation matrix with osci. freq
            rho = .9;
            theta = (2*pi*osc_freq)*(1/fs);
            for i=1:k
                for j=1:k
                    if i == j
                        a0(i*2-1,j*2-1) = rho*cos(theta);
                        a0(i*2,j*2-1) = rho*sin(theta);
                        a0(i*2-1,j*2) = rho * -sin(theta);
                        a0(i*2,j*2) = rho * cos(theta);                        
                    else
                        rdrho = a + (b-a)*rand(1);
                        rdtheta = -pi + (2*pi)*rand(1);
                        
                        a0(i*2-1,j*2-1) = rdrho * cos(rdtheta);
                        a0(i*2,j*2-1) = rdrho * sin(rdtheta);
                        a0(i*2-1,j*2) = rdrho * -sin(rdtheta);
                        a0(i*2,j*2) = rdrho * cos(rdtheta);
                    end
                end
            end
        end
        
        % ===== create random A matrix =====
        % randomly weakly correlated
        % rho between 0.01~0.2 -- a=0.01 & b=0.2
        % random phase between -pi to pi
        % random value between interval (a,b) using a + (b-a)*rand(1)
        function a0 = randa2(k,osc_freq,fs,aAR)
            a0 = eye(2*k);
            a = 0.01;
            b = 0.2;
            % on diagonal - rotation matrix with osci. freq
            rho = aAR;
            theta = (2*pi*osc_freq)*(1/fs);
            for i=1:k
                for j=1:k
                    if i == j
                        a0(i*2-1,j*2-1) = rho*cos(theta);
                        a0(i*2,j*2-1) = rho*sin(theta);
                        a0(i*2-1,j*2) = rho * -sin(theta);
                        a0(i*2,j*2) = rho * cos(theta);                        
                    else
                        rdrho = a + (b-a)*rand(1);
                        rdtheta = -pi + (2*pi)*rand(1);
                        
                        a0(i*2-1,j*2-1) = rdrho * cos(rdtheta);
                        a0(i*2,j*2-1) = rdrho * sin(rdtheta);
                        a0(i*2-1,j*2) = rdrho * -sin(rdtheta);
                        a0(i*2,j*2) = rdrho * cos(rdtheta);
                    end
                end
            end
        end
        
        % ===== visualize Q matrix - magnitude & phase =====
        function [mag, phase] = Q_mag_phase(q)
            n = size(q,1)/2;
            mag = zeros(n,n);
            phase = nan(n);
            
            for i=1:n
                for j=1:n
                    if i ~= j
                        mag(i,j) = sqrt(q(i*2-1,j*2-1)^2+q(i*2-1,j*2)^2);
                    end
                end
            end
            
            for i=1:n
                for j=1:n
                    if i ~= j
                        if mag(i,j) > 0.15
                            phase(i,j) = atan2(q(i*2,j*2-1),q(i*2,j*2));
                        end
                    end
                end
            end
        end
        
        % ===== visualize B matrix - magnitude & phase =====
        function [mag, phase] = B_mag_phase(b)
            n = size(b,1);
            k = size(b,2)/2;
            
            mag = zeros(n,k);
            phase = nan(n,n,k);
            
            % magnitude
            for i=1:n %node
                for j=1:k %oscillator
                    mag(i,j) = sqrt(b(i,j*2-1)^2+b(i,j*2)^2);
                end
            end
            % phases
            for g=1:k
                for i=1:n
                    for j=1:n
                        if i==j
                            phase(i,j,g) = NaN;
                        else
                            v1 = [b(i,g*2-1), b(i,g*2)];
                            v2 = [b(j,g*2-1), b(j,g*2)];
                            ag1 = atan2(v1(1),v1(2));
                            ag2 = atan2(v2(1),v2(2));
                            if (mag(i,g)>0.174 && mag(j,g)>0.174)
                                if atan2(v1(1),v1(2))-atan2(v2(1),v2(2))>=0
                                    phase(i,j,g) = -mod(atan2(v1(1),v1(2))-atan2(v2(1),v2(2)),2*pi);
                                else
                                    phase(i,j,g) = mod((atan2(v2(1),v2(2))-atan2(v1(1),v1(2))),2*pi);
                                end
                            end
                        end
                    end
                end
            end
        end

        % ===== visualize A matrix - magnitude & phase =====
        function [mag, phase] = A_mag_phase(a)
            n = size(a,1)/2;
            mag = zeros(n,n);
            phase = zeros(n,n);
            threshold = 0.15;
            for i=1:n
                for j=1:n
                    if i == j
                        mag(i,j) = 0;    
                    else
                        mag(i,j) = sqrt(a(i*2-1,j*2-1)^2+a(i*2-1,j*2)^2);
                    end
                end
            end
            
            for i=1:n
                for j=1:n
                    if i == j
                        phase(i,j) = NaN;     
                    else
                        if mag(i,j) < threshold
                            phase(i,j) = NaN;
                        else
                            phase(i,j) = atan2(a(i*2,j*2-1),a(i*2,j*2));
                        end
                    end
                end
            end
        end

        % ===== handle plotting phase =====
        function h = imagesc2(img)
            h = imagesc(img);
            axis image off
            if ndims(img) == 2
                set(h, 'AlphaData', ~isnan(img))
            elseif ndims( img ) == 3
                set(h, 'AlphaData', ~isnan(img(:, :, 1)))
            end
            if nargout < 1
                clear h
            end
        end
    end
end

    