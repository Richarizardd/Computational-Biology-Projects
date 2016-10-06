function simulateSystem(x0, A, B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, goal, num1, num2, TITLE, K, G)
    %%% Simulation
    % Preallocation.
    xhat = cell(1,pSteps);
    x = cell(1,pSteps);
    y = cell(1, pSteps);
    u = cell(1, pSteps);
    
    for k = 1:pSteps
        xhat{k} = zeros(50,7);
        x{k} = zeros(50,7);
    end
    
    for ind = 1:numAverage
        % Initialization.
        xhat{1}(ind,:) = x0;
        x{1}(ind,:) = x0;
    end
    
    y{1} = H*x{1}(ind,:)'+mvnrnd(zeros(2,1),Qy)';
    u{1} = -G{1}*xhat{1}(ind,:)';
    
    for ind = 1:numAverage
        for k = 1:pSteps-1
            x{k+1}(ind,:) = A*x{k}(ind,:)'+B*u{k}+...
                mvnrnd(zeros(7, 1), Qx)'+B*(C1*u{k}*randn+C2*u{k}*randn);

            if k < Thold
                x{k+1}(ind,4:5) = 0;
            end

            y{k+1} = H*x{k}(ind,:)' + mvnrnd(zeros(2, 1), Qy)';
            xhat{k+1}(ind,:) = A*xhat{k}(ind,:)'+...
                A*K{k}*(y{k}-H*xhat{k}(ind,:)')+B*u{k};
            u{k+1} = -G{k+1}*xhat{k+1}(ind,:)';
        end
    end
    
    eye_pos = zeros(1,pSteps);
    head_pos = zeros(1,pSteps);
    gaze_pos = zeros(1,pSteps);
    
    for k = 1:pSteps
        temp_avg = x{k}(1,:);
        eye_pos(k) = temp_avg(1);
        head_pos(k) = temp_avg(4);
        gaze_pos(k) = temp_avg(1)+temp_avg(4);
    end
    
    figure(num1);
    plot(Delta*(1:pSteps), eye_pos)
    hold on
    plot(Delta*(1:pSteps), head_pos)
    plot(Delta*(1:pSteps), gaze_pos)
    plot(get(gca,'xlim'), [goal goal])
    xlabel({'Seconds (sec)'},'FontSize',12);
    ylabel({'Radians'},'FontSize',12);
    title({TITLE},'FontSize',16);
    legend('Eye','Head','Gaze', 'Goal')

    for k = 1:pSteps
        temp_avg = mean(x{k});
        eye_pos(k) = temp_avg(1);
        head_pos(k) = temp_avg(4);
        gaze_pos(k) = temp_avg(1)+temp_avg(4);
    end
    
    figure(num2);
    plot(Delta*(1:pSteps), eye_pos)
    hold on
    plot(Delta*(1:pSteps), head_pos)
    plot(Delta*(1:pSteps), gaze_pos)
    plot(get(gca,'xlim'), [goal goal])
    xlabel({'Seconds (sec)'},'FontSize',12);
    ylabel({'Radians'},'FontSize',12);
    title([TITLE ', Averaged'],'FontSize',16);
    legend('Eye','Head','Gaze', 'Goal')
end