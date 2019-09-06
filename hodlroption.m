function opt = hodlroption(key, value)
%HMOPTION Set or get an option for the hodlr toolbox.
%
% Valid options are:
%   'block-size': Integer representing the minimum block size.
%   'threshold': Value used for off-diagonal truncation.
%   'dense-compression': Can be either 'rrqr' or 'svd', and selects the
%       method used to compress unstructured dense matrices when calling
%       hodlrA = hodlr(A);

global hodlr_block_size
global hodlr_threshold
global hodlr_dense_compression

if isempty(hodlr_dense_compression)
	hodlr_dense_compression = 'rrqr';
end

if isempty(hodlr_block_size)
    hodlr_block_size = 256;
end

if isempty(hodlr_threshold)
    hodlr_threshold = 1e-12;
end

if ~exist('key', 'var')
    error('Please specify a key');
end

if ~exist('value', 'var')
    switch key
        case 'block-size'
            opt = hodlr_block_size;
        case 'threshold'
            opt = hodlr_threshold;
		case 'dense-compression'
			opt = hodlr_dense_compression;
        otherwise
            error('Unsupported option specified');
    end
else
    switch key
        case 'block-size'
            if value <= 2
                error('minimum block size must be at least 3');
            else
                hodlr_block_size = value;
            end
        case 'threshold'
            if value < 0
                error('threshold has to be positive');
            else
                hodlr_threshold = max(eps, value);
			end
		case 'dense-compression'
			if ~strcmp(value, 'rrqr') && ~strcmp(value, 'svd') && ...
					~strcmp(value, 'lanczos')
				error('Invalid value for dense-compression');
			else
				hodlr_dense_compression = value;
			end
        otherwise
            error('Unsupported option specified');
    end
    
end
