function out = softEq(A, B, precision, dim)

%
%   SOFT EQUALS
%
%  EXAMPLE: softEquals = softEq(A, B, precision, dim);
%
% Determines the element-wise equality between matricies A and B but
% neglects any round off error. The amount error to be tolerated is
% specified by the optional 'precision' argument. The default tolerance is
% ten decimal places. The optional argument `dim` specifies whether the
% equality should be computed along the 'rows' or 'cols'.
%
% CAH 03/09

if nargin < 2 || isempty(A) || isempty(B)
    out = [];
    return
end

if nargin < 3 || isempty(precision);
    precision = 10;
end

%if the optional argument dim is specified then treat this as a call to
%ismember
if nargin > 3
    if strcmpi(dim, 'cols')
        A = A';
        B = B';
    end
    
    %make sure A is the smaller of the two matricies. I'll use A as the
    %template and look into B for instances of A.
    if size(A,1) > 1
        error(['the template (first argument) must be a row or column vector' ...
            ' (depending on what you passed in as `dim`)'])
    end
    matches = abs(B - repmat(A, size(B, 1), 1)) < 10^(-precision);
    out = sum(matches, 2) == size(B, 2);
    
    %transform back if the caller specified columns
    if strcmpi(dim, 'cols')
        out = out';
    end
else
    %just the element wise soft equals:
    out = abs(A-B) < 10^(-precision);
end
