function[A, B] = calculateAandB(tau1eye, tau2eye, alpha1eye, tau1head, tau2head, alpha1head , Delta)
    alpha2 = 1;
    
    % parameters of the eye
    k_e = 1;
    b_e = tau1eye+tau2eye;
    m_e = tau1eye*tau2eye;

    % parameters of the head
    b_h = tau1head+tau2head;
    m_h = tau1head*tau2head;
    
    % third-order system
    A_e = [0 1 0; -k_e/m_e -b_e/m_e 1/m_e; 0 0 -alpha2/alpha1eye];
    B_e = [0; 0; 1/alpha1eye];
    A_h = [0 1 0; -k_e/m_h -b_h/m_h 1/m_h; 0 0 -alpha2/alpha1head];
    B_h = [0; 0; 1/alpha1head];
    
    % noise-free version of dynamical system in continuous time
    A_c = vertcat(horzcat(A_e, zeros(3)), horzcat(zeros(3), A_h));
    B_c = vertcat(horzcat(B_e, zeros(3, 1)), horzcat(zeros(3, 1), B_h));
    
    A = vertcat(horzcat(expm(A_c*Delta), zeros(6, 1)), horzcat(zeros(1, 6), 1));
    B = vertcat(inv(A_c)*(expm(A_c*Delta)-eye(6))*B_c, zeros(1, 2));
end