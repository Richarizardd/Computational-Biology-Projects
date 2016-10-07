%% HW10 - 1. Computer experiment with a Hopfield network
% In this experiment, I built a Hopfield network with 1600 binary neurons
% that are fully connected with one another, and examined the overall
% behaviors of a current network.

%% 0 - Initialization
load HOP
W = A*A' + B*B' + C*C';

figure('uni','pi','pos',[100 100 1200 300]);
Patterns = [A B C D];
Names = ['A' 'B' 'C' 'D'];

for ind = 1:4
    subplot(1,4,ind)
    imagesc(reshape(Patterns(:,ind),40,40))
    colormap(gray)
    title(['Pattern ' Names(ind)])
    axis('square');
end

%% 1 - Stored Memories
% All of these states are stationary, as the state S (representing A, B, or 
% C) was equivalent to sign(W*S) after one iteration. This makes sense,
% because form using the weight rule (W = A*A' + B*B' + C*C'), each memory
% pattern becomes an eigenvector of the weight matrix, and thus, A, B, and
% C are already memory patterns in our network. By starting at a state
% identical to one stored in memory, our network will instantly converge to
% that memory.

%%
disp(isequal(A, sign(W*A)));
disp(isequal(B, sign(W*B)));
disp(isequal(C, sign(W*C)));

figure('uni','pi','pos',[100 100 1200 300]);
Patterns = [A B C];
Names = ['A' 'B' 'C'];

for ind = 1:3
    S = Patterns(:,ind);
    subplot(1,3,ind)
    imagesc(reshape(sign(W*S),40,40))
    colormap(gray)
    title(['State ' Names(ind)])
    axis('square');
end

%% 2 - Noise in Initial State
% Here, I observed that Test1, Test2, and Test3 converged to the same final
% stationary state (Pattern C). This makes sense, as the pattern of Test1
% is a rotation of Pattern C. Since our network starts from a state
% sufficiently close to one already stored in memory, it will reach that
% memory in the next time step. Test2 and Test3 is just Test1 that have had
% random weights flipped, so as a result, Pattern C is still a local
% minimum in the Liaponov function of our Hopfield.


%%
figure('uni','pi','pos',[100 100 900 300]);
Tests = [Test1 Test2 Test3];

for ind = 1:3
    S = Tests(:,ind);
    stationary = 0;
    
    while stationary == 0
        tempS = sign(W*S);
        if isequal(S, tempS)
            stationary = 1;
        end
        S = tempS;
    end
    
    subplot(1,3,ind);
    imagesc(reshape(S,40,40))
    colormap(gray)
    title(['Final Stationary State of Test' int2str(ind)])
    axis('square');
end


%% 3 - Unrelated pattern
% Final state for pattern D. The synaptic weight network was unable to
% remember the lochness pattern in D, as it converged to Pattern B. This is
% because Pattern D was never stored as a memory in our network, and when
% we use it in our input, our network will try and converge it to the
% closest memory state, which happened to be Pattern B.

%%
S = D;
stationary = 0;

while stationary == 0
    tempS = sign(W*S);
    if isequal(S, tempS)
        stationary = 1;
    end
    S = tempS;
end

figure
imagesc(reshape(S,40,40))
colormap(gray)
title('Final Stationary State of D after before D to W')

%% 4 - Learning new memory pattern
% After adding D into the synpatic weight network, we were able to recover
% the lochness pattern (Pattern D). This makes sense since we added Pattern
% D as a memory state in our network.

%%
W = W + D*D';
S = D;
stationary = 0;

while stationary == 0
    tempS = sign(W*S);
    if isequal(S, tempS)
        stationary = 1;
    end
    S = tempS;
end

figure
imagesc(reshape(S,40,40))
colormap(gray)
title('Final Stationary State of D after adding D to W')

%% 5 - Are old memories affected?
% Old memories were unaffected, as Test3 converged to the final stationary
% state found in Question 2 (Pattern C). Pattern C is still a local minimum
% in the Liapnov function of our Hopfield, and different memory patterns
% that are statistically independent of one another can be stored in the
% same network without interferences.

%%
S = Test3;
stationary = 0;
while stationary == 0
    tempS = sign(W*S);
    if isequal(S, tempS)
        stationary = 1;
    end
    S = tempS;
end

figure
imagesc(reshape(S,40,40))
colormap(gray)
title('Final Stationary State of Test3 after adding D to W')

%% 6 - Spurious memories
% In setting random initial states, I found that not all of the states
% converged to the same final stationary state. Different seeds will
% converge to states in our network, some of which are not real memory
% states due to interference if an initial state is statistically close to
% multiple memories.

%%
figure('uni','pi','pos',[100 100 900 900]);

for i = 1:16
    S = sign(randn(1600,1));
    stationary = 0;
    
    while stationary == 0
        tempS = sign(W*S);
        if isequal(S, tempS)
            stationary = 1;
        end
        S = tempS;
    end
    
    subplot(4,4,i)
    imagesc(reshape(S,40,40))
    colormap(gray)
end

%% 7 - Symmetric connections
% By multipying each term (Patterns A, B, C, and D) with its tranpose, each
% resulting matrix will be symmetric and square. As a result, the weight
% rule does force the weight matrix to be symmetric as well. A more formal
% proof is given below.
%%
% 
% <<10_1_7_Proof.png>>
% 

%%
figure('uni','pi','pos',[100 100 1000 500]);
subplot(1,2,1)
imagesc(W)
colormap(jet)
title('W')

subplot(1,2,2)
imagesc(W(1:100, 1:100))
colormap(jet)
title('W(1:100,1:100)')

%% 8 - Damage to synaptic connections
% Though a majority of connections were destroyed, Test3 still converged to
% the same final stationary network as found in Question 2 (Pattern C).
% This is due to the robust nature of Hopfields, as contribution from any 
% damaged, but also small subsets of neurons is relatively small, and thus 
% can be completely ignored without affecting the state of the network in 
% the next time step.

%%
W=randzero(W, 0.7);
S = Test3;
stationary = 0;

while stationary == 0
    tempS = sign(W*S);
    if isequal(S, tempS)
        stationary = 1;
    end
    S = tempS;
end

figure
imagesc(reshape(S,40,40))
colormap(gray)
title('Final Stationary State of Test3 after destorying connections in W')