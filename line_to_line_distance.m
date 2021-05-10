function [d] = line_to_line_distance(M1, u1, M2, u2, verbose)
%% line_to_line_distance : function to compute the minimum distance between
% two lines L1(M1,u1), and L2(M2,u2) of the 3D space.
%
% Author & support nicolas.douillet (at) free.fr, 2020.
%
%
% Syntax
%
% d = line_to_line_distance(M1, u1, M2, u2);
% d = line_to_line_distance(M1, u1, M2, u2, verbose);
%
%
% Description
%
% d = line_to_line_distance(M1, u1, M2, u2) computes the distance between
% L1 and L2 and display a message in case they are parallel.
%
% d = line_to_line_distance(M1, u1, M2, u2, verbose) computes the distance between
% L1 and L2 and display a message in case they are parallel if verbose
% is set to logical *true/1.  
%
%
% Input arguments
%
% - M1 : real row or column vector double, a point belonging to L1, numel(M1) = 3.
% - u1 : real row or column vector double, one L1 director vector, numel(u1) = 3.
% - M2 : real row or column vector double, a point belonging to L2, numel(M2) = 3.
% - u2 : real row or column vector double, one L2 director vector, numel(u2) = 3.
% - verbose : either logical *true/false or numeric *1/0.
%
%
% Output argument
%
% - d : real scalar double, the -minimum- distance between the two lines.
%
%
% Example #1
%
% u1 = [0 2 -1];
% M1 = [1 2 -1];
% u2 = [0 2 1];
% M2 = [-1 2 1];
% d = line_to_line_distance(M1, u1, M2, u2) % d = 2 expected 
%
%
% Example #2
%
% u1 = [-1 1 0]';
% M1 = [1 1 1]';
% u2 = [1 1 -2]';
% M2 = [-1 -1 -1]';
% d = line_to_line_distance(M1, u1, M2, u2) % d = 2*sqrt(3) expected 
%
%
% Example #3 : parallel lines
%
% u1 = [-1 1 0];
% M1 = ones(1,3)/3;
% u2 = [1 -1 0]; % = -u1
% M2 = -2*ones(1,3)/3;
% d = line_to_line_distance(M1, u1, M2, u2) % d = sqrt(3) expected, L1 // L2
%
%
% Example #4 : intersecting lines
%
% u1 = [0 1 0]';
% M1 = [0 2 0]';
% u2 = [0 4 -7]';
% M2 = u2;
% d = line_to_line_distance(M1, u1, M2, u2) % d = 0 expected
%% Input parsing
assert(nargin > 3,'Not enough input argument.');
assert(nargin < 6,'Too many input arguments.');
if nargin < 5   
    verbose = true;    
end
assert(isequal(size(u1),size(u2),size(M1),size(M2)),'All inputs vectors and points must have the same size.');
assert(isequal(numel(u1),numel(u2),numel(M1),numel(M2),3),'All inputs vectors and points must have the same number of elements (3).');
assert(isequal(ndims(u1),ndims(u2),ndims(M1),ndims(M2),2),'All inputs vectors and points must have the same number of dimensions (2).');
assert(isreal(u1) && isreal(u2) && isreal(M1) && isreal(M2),'All inputs vectors and points must contain real numbers only.');
%% Body
epsilon = eps;
M1M2 = M2 - M1;
w = cross(u1,u2);
nw = norm(w);
if norm(w) < 1e-4*epsilon % case L1 // L2
    if verbose
       
        disp('Lines (M1,u1) and (M2,u2) are parallel.');
        
    end
    
    d_M1H = dot(M1M2,u1/norm(u1));
    
    d = sqrt(norm(M1M2)^2 - d_M1H^2); % = dM2H (Pythagore's theorem)
    
else
        
    d = dot(w,M1M2)/nw;
        
end
if verbose
    
    if d == 0
        
        disp('Lines (M1,u1) and (M2,u2) intersect in at least one point.');
        
        if nw < 1e-4*epsilon
            
            disp('Lines (M1,u1) and (M2,u2) also have colinear direct vectors : they are one same line.');
            
        end
        
    end
    
end
end % line_to_line_distance