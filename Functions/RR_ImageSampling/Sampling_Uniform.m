function [ varargout ] = Sampling_Uniform( image_, nSamples_, displaySamples_ )
%SAMPLING_UNIFORM Spatially samples an image using a uniform distribution.
%   INPUTS:
%       image_ - image that will be sampled (M x N x D matrix)
%       nSamples_ - number of desired samples (default: 50)
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
    validateattributes( nSamples_, {'numeric'}, {'positive','integer','numel',1} );
end

if nargin < 3
    displaySamples_ = false;
else
    validateattributes( displaySamples_, {'logical'}, {'numel',1} );
end

%% generate sample points
% first row are the samples along the columns, second row are the samples
% down the rows
samplePoints = [ randi( size( image_, 2 ), [ 1 nSamples_ ] );
                    randi( size( image_, 1 ), [ 1 nSamples_ ] ) ];

%% make sure points lie within the image boundary
samplePoints(1,:) = min( max( 1, samplePoints( 1, : ) ), size( image_, 2 ) );
samplePoints(2,:) = min( max( 1, samplePoints( 2, : ) ), size( image_, 1 ) );
samplePoints = ceil( samplePoints );
varargout{1} = samplePoints;

%% obtain the values from the image
if nargout > 1
    sampleValues = zeros( size( image_, 3 ), nSamples_ );
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
