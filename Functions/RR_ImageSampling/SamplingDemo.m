%% demonstration of the sampling algorithms
imageExample = zeros( 640, 480, 50 ); % example, multidimensional image

% run the sampling algorithms sampling 200 points
[ p_g, v_g ] = Sampling_Grid( imageExample, 200, true );
[ p_j, v_j ] = Sampling_Jittered( imageExample, 200, true );
[ p_u, v_u ] = Sampling_Uniform( imageExample, 200, true );
% use 50 candidates for the best candidate algorithm
[ p_bc, v_bc ] = Sampling_BestCandidate( imageExample, 200, 50, true );