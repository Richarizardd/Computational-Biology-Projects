function G = calculateFeedbackGains(x0, A, B, pSteps, k1Steps, Qx, Qy, ...
    H, C1, C2, K, T, L)
    
    % init
    W_x{pSteps} = H'*T*H;
    W_e{pSteps} = zeros(7);
    C_x{pSteps} = C1'*B'*W_x{pSteps}*B*C1 + C2'*B'*W_x{pSteps}*B*C2;
    C_e{pSteps} = C1'*B'*W_e{pSteps}*B*C1 + C2'*B'*W_e{pSteps}*B*C2;
    G{pSteps} = zeros(2,7);
    G{pSteps-1} = inv(L+C_x{pSteps}+...
        C_e{pSteps}+B'*W_x{pSteps}*B)*B'*W_x{pSteps}*A;
    
    for k = pSteps-1:-1:2
        
        if k <= k1Steps
            T = zeros(2);
        end
        
        % Update parameters for optimal policy
        W_e{k} = (A-A*K{k}*H)'*W_e{k+1}*(A-A*K{k}*H) + ...
            G{k}'*B'*W_x{k+1}*A;
        W_x{k} = H'*T*H + A'*W_x{k+1}*A - G{k}'*B'*W_x{k+1}*A; % D = 0
        
        % Update parameters for optimal motor command
        C_x{k} = C1'*B'*W_x{k}*B*C1 + C2'*B'*W_x{k}*B*C2;
        C_e{k} = C1'*B'*W_e{k}*B*C1 + C2'*B'*W_e{k}*B*C2;
        
        % Derive optimal control opolicy at time step k-1
        G{k-1} = inv(L+C_x{k}+C_e{k}+B'*W_x{k}*B)*B'*W_x{k}*A;
    end
end