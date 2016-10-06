function K = calculateKalmanGains(x0, A, B, pSteps, Qx, Qy, H, C1, C2, G)
    S_e{1} = zeros(7); % (1|0) - posterior
    S_x{1} = x0*x0'; % (1|0) - posterior
    S_xe{1} = zeros(7); % (1|0) - posterior
    K{1} = S_e{1}*H'*inv(H*S_e{1}*H'+Qy); % 
    
    for k = 2:pSteps
        S_e{k} = A*(eye(7)-K{k-1}*H)*S_e{k-1}*A'+Qx+...
            B*C1*G{k-1}*S_x{k-1}*G{k-1}'*C1'*B'+...
            B*C2*G{k-1}*S_x{k-1}*G{k-1}'*C2'*B';
        S_x{k} = (A-B*G{k-1})*S_x{k-1}*(A-B*G{k-1})'+...
            (A-B*G{k-1})*S_xe{k-1}*(A*K{k-1}*H)'+...
            A*K{k-1}*H*S_xe{k-1}'*(A-B*G{k-1})'+...
            A*K{k-1}*H*S_e{k-1}*A';
        S_xe{k} = (A-B*G{k-1})*S_xe{k-1}*(eye(7)-K{k-1}*H)'*A';
        K{k} = S_e{k}*H'*inv(H*S_e{k}*H'+Qy);
    end  
end

