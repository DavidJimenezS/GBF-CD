function [ varargout ] = Sampling_BestCandidate( image_, nSamples_, nCandidates_, displaySamples_ )
%SAMPLING_BESTCANDIDATE Sample points using Mitchell's best candidate
%sampling. Slower than other methods but provides greater spatial coverage
%   INPUTS:
%       image_ - image that will be sampled (M x N x D matrix)
%       nSamples_ - number of desired samples (default: 50). 
%       nCandidates_ - number of candidates generated for each sample point
%       (default: 30)
%       displaySamples_ - show the samples obtained on the image (default:
%       false)
%
%   OUTPUTS:
%       samplePoints - positions at which the samples were obtained ([x;y])
%       sampleValues - (optional) pixel intensites at the samplePoints (D x nSamples_)
%
%   @author Rishi Ramakrishnan
%   @version 1.0
%   @date 5 December 2015

%% validate the input paramters
if nargin < 2 % user has not inputted the number of samples so set the default value
    nSamples_ = 50;
else
    if isempty( nSamples_ )
        nSamples_ = 50;
    end
    validateattributes( nSamples_, {'numeric'}, {'positive','integer','numel',1} );
end

if nargin < 3
    nCandidates_ = 30;
else
    if isempty( nCandidates_ )
        nCandidates_ = 30;
    end
    validateattributes( nCandidates_, {'numeric'}, {'positive','integer','numel',1} );
end

if nargin < 4
    displaySamples_ = false;
else
    if isempty( displaySamples_ )
        displaySamples_ = false;
    end
    validateattributes( displaySamples_, {'logical'}, {'numel',1} );
end

%% generate sample points
samplePoints = [ randi( size( image_, 2 ), 1 ); randi( size( image_, 1 ), 1 ) ];
for i = 1:1:nSamples_
    candidates = [ rand( [ 1 nCandidates_ ] )*size(image_,2); rand( [ 1 nCandidates_ ] )*size(image_,1) ];
    D = pdist2( samplePoints', candidates' ); % calculate distance between candidates and point in the sample set
    D = min( D, [], 1 ); % get the minimum distance to the points already in the set
    [ Y, I ] = max( D ); % get the maximum of the distances and add to the set
    samplePoints( :, end+1 ) = ceil( candidates( :, I(1) ) );
end
% make sure the points are within the image boundaries
samplePoints(1,:) = min( max( 1, samplePoints( 1, : ) ), size( image_, 2 ) );
samplePoints(2,:) = min( max( 1, samplePoints( 2, : ) ), size( image_, 1 ) );

varargout{1} = samplePoints;
%% obtain the values from the image
if nargout > 1
    sampleValues = zeros( size( image_, 3 ), size( samplePoints, 2 ) );
    for i = 1:1:size( image_, 3 )
        temp = image_( :, :, i );
        sampleValues( i, : ) = temp( sub2ind( size( temp ), samplePoints( 2, : ), samplePoints( 1, : ) ) );
    end
    varargout{2} = sampleValues;
end
%% display the samples if required
if displaySamples_
    figure
    imagesc( image_( :, :, 1 ) )
    hold on
    plot( samplePoints( 1, : ), samplePoints( 2, : ), 'rx', 'LineWidth', 2 )
end


end
