function[K,G] = simulateCase(x0, A, B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, goal, num1, num2, TITLE)
    
    G = cell(1, pSteps);
    for ind = 1:pSteps
        G{ind} = zeros(2,7);
    end
    
    
    for ind = 1:10
        K = calculateKalmanGains(x0, A, B, pSteps, Qx, Qy, H, C1, C2, G);
        G = calculateFeedbackGains(x0, A, B, pSteps, k1Steps, Qx, Qy, ...
            H, C1, C2, K, T, L);
    end
    
    simulateSystem(x0, A, B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, goal, num1, num2, TITLE, K, G);
end