function [ varargout ] = Sampling_Stratified( image_, nSamples_, displaySamples_ )
%SAMPLING_STRATIFIED Spatially samples an image using a stratified approach.
%   INPUTS:
%       image_ - image that will be sampled (M x N x D matrix)
%       nSamples_ - number of desired samples (default: 50). Can be a 1
%       [ samplesPerImage ] or 2 element [ samplesPerRow, samplesPerColumn ] vector.
%       displaySamples_ - show the samples obtained on the image (default:
%       false)
%
%   OUTPUTS:
%       samplePoints - positions at which the samples were obtained ([x;y])
%       sampleValues - (optional) pixel intensites at the samplePoints (D x nSamples_)
%
%   ALGORITHM:
%       Splits the image space into strata and then samples each strata
%       using a uniform distribution.
%
%   @author Rishi Ramakrishnan
%   @version 2.0
%   @date 6 December 2015

%% validate the input paramters
if nargin < 2 % user has not inputted the number of samples so set the default value
    nSamples_ = 50;
else
    if isempty( nSamples_ )
        nSamples_ = 50;
    end
    validateattributes( nSamples_, {'numeric'}, {'positive','integer'} );
end
if numel( nSamples_ ) == 1 
    % calculate how many row and column samples are required to generate a total of approximately nSamples_ samples
    % number of samples for the rows and columns are in the same ratio as
    % the image
    nSamplesCol = size( image_, 2 ) * sqrt( nSamples_ / numel( image_( :, :, 1 ) ) );
    nSamplesRow = round( size( image_, 1 ) / size( image_, 2 ) * nSamplesCol );
    nSamplesCol = round( nSamplesCol );
elseif numel( nSamples_ ) == 2
    [ nSamplesRow, nSamplesCol ] = deal( nSamples_(1), nSamples_(2) );
else
    display( 'Sampling_Stratified: nSamples_ must be a vector of 1 or 2 elements' )
    return;
end

if nargin < 3
    displaySamples_ = false;
else
    if isempty( displaySamples_ )
        displaySamples_ = false;
    end
    validateattributes( displaySamples_, {'logical'}, {'numel',1} );
end

%% generate sample points
% samples are generated with regular spacing by calculating the even
% spacing between each point along the columns and down the row
dCol = ( size( image_, 2 ) / ( nSamplesCol + 1 ) );
dRow = ( size( image_, 1 ) / ( nSamplesRow + 1 ) );

% start of each strata
samplePointsCol = 0:dCol:( size( image_, 2 ) - dCol );
samplePointsRow = 0:dRow:( size( image_, 1 ) - dRow);

% generate all the sample points
nSamplesCol = numel( samplePointsCol );
nSamplesRow = numel( samplePointsRow );
samplePointsRow = repmat( samplePointsRow', 1, nSamplesCol );
samplePointsCol = repmat( samplePointsCol, nSamplesRow, 1 );

samplePoints = [ reshape( samplePointsCol, 1, numel( samplePointsCol ) ) + rand( 1, numel( samplePointsCol ) )*dCol;
                    reshape( samplePointsRow, 1, numel( samplePointsRow ) ) + rand( 1, numel( samplePointsRow ) )*dRow ];

%% make sure points lie within the image boundary
samplePoints(1,:) = min( max( 1, samplePoints( 1, : ) ), size( image_, 2 ) );
samplePoints(2,:) = min( max( 1, samplePoints( 2, : ) ), size( image_, 1 ) );
samplePoints = ceil( samplePoints );
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
