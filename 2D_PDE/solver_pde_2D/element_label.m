% ---------------------------------------------------- %
% function [elab] = element_label(lx,ly)
%
% inputs: lx = level of refinement in x direction
%         ly = level of refinement in y direction
%
% output: elab = element label for every grid points
% ---------------------------------------------------- %
function [elab] = element_label(lx,ly)
mel=0; % macro-element consist of 4 elements
mx = 2^(lx-1); % number of macro-element in x direction
my = 2^(ly-1); % number of macro-element in y direction
nx = 2^lx;     % number of element in x direction
kx = 1;
ky = 1;
for j=1:my
   for i=1:mx
      mref=(nx+1)*(ky-1)+kx;
      mel=mel+1;
      nvv(1) = mref;
      nvv(2) = mref+2;
      nvv(3) = mref+2*nx+4;
      nvv(4) = mref+2*nx+2;
      nvv(5) = mref+1;
      nvv(6) = mref+nx+3; 
      nvv(7) = mref+2*nx+3; 
      nvv(8) = mref+nx+1;
      nvv(9) = mref+nx+2; 
      mv(mel,1:9)=nvv(1:9);
      kx = kx + 2;
   end
   ky = ky + 2; 
   kx = 1;
end


for k=1:mel
% first element
ke = 4*k-3;
elab(ke,1) = mv(k,1);
elab(ke,2) = mv(k,5);
elab(ke,3) = mv(k,9);
elab(ke,4) = mv(k,8);
% second element
ke = 4*k-2;
elab(ke,1) = mv(k,5);
elab(ke,2) = mv(k,2);
elab(ke,3) = mv(k,6);
elab(ke,4) = mv(k,9);
% third element
ke = 4*k-1;
elab(ke,1) = mv(k,9);
elab(ke,2) = mv(k,6);
elab(ke,3) = mv(k,3);
elab(ke,4) = mv(k,7);
% fourth element
ke = 4*k;
elab(ke,1) = mv(k,8);
elab(ke,2) = mv(k,9);
elab(ke,3) = mv(k,7);
elab(ke,4) = mv(k,4);
end

end

