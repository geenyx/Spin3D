function [R] = eul2rot(seq, varargin)
% seq is the order / sequence of rotations, e.g, 'xyz', 'xz', 'yz'. 
% For dual axis tomography, example input would be:
% eul2rot('zy',[0,30],1:180)
% which reads: z tilts are 0 and 30. After a tilt is applied it is followed
% by a y rotation. Projetions are taken for each tilt and y combination.
% The y rotation varies between 1 and 180 degrees
% varargin are vectors e.g, [0:90], [1,2]
% the output is the combination of the size of all angles N1 x N2 x N3 etc
%   this script probably doesnt scale well with large inputs - not
%   vectorised

if size(varargin,2) ~= size(seq,2)
    error("Number of axes do not match number of rotations")
end

D = varargin;
[D{:}] = ndgrid(varargin{:});
vector = cell2mat(cellfun(@(m)m(:),D,'uni',0));

nseq = size(vector,2);
nangles = size(vector,1);

% initialize R to identity
R = repmat(eye(3,3),[1,1,nangles]);

for ii=1:nseq
    s = seq(ii);
    switch s
        case 'x'
            for jj=1:nangles
                R(:,:,jj) = rotx(vector(jj,ii))*R(:,:,jj);
            end
        case 'y'
            for jj=1:nangles
                R(:,:,jj) = roty(vector(jj,ii))*R(:,:,jj);
            end
        case 'z'
            for jj=1:nangles
                R(:,:,jj) = rotz(vector(jj,ii))*R(:,:,jj);
            end
        otherwise
            error("Axis not defined")
    end
end
