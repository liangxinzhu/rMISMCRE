% ---------------------------------------------------------------- %
% function [varargout] = LGCprior(lx,ly,rate,a,b,mu,type)
%
% inputs: lx   = level of refinement in x direction
%         ly   = level of refinement in y direction
%         rate = \beta^{\prime}
%         a    = \theta_2
%         b    = \theta_3
%         mu   = \theta_1
%         type = 1: single level
%         type = 2: increments of increments along x
%         type = 3: increments of increments along y
%         type = 4: increments of increments
%         type = 5: increments
%
% outputs: varargout{1} = u1
%          varargout{1} = u2
%          varargout{1} = u3
%          varargout{1} = u4
% ---------------------------------------------------------------- %

function [varargout] = LGCprior(lx,ly,rate,a,b,mu,type)
Nx = 2^lx; Ny = 2^ly;
kx = -Nx/2:Nx/2-1; ky = -Ny/2:Ny/2-1;
[k1,k2]=meshgrid(kx,ky);
% complex white noise in Fourier domain
frr=fftshift(fft2(randn(Ny,Nx)))/sqrt(Nx*Ny);
% sqrtcov1 = sqrt(a)*sqrt(k1.^2+k2.^2+b).^(-rate/2);
% sqrt covariance
sqrtcov1 = sqrt(a)*sqrt(k1.^2+b).^(-rate/2).* ...
        sqrt(k2.^2+b).^(-rate/2);
switch type
    case 1 % single level
        fxi1 = sqrtcov1.*frr;
        % shift back to spatial domain
        samp1 = real(ifft2(ifftshift(fxi1)))*Nx*Ny;
        varargout{1} = mu + samp1;
    case 2 % increments of increments along x
        fxi1 = sqrtcov1.*frr;
        samp1 = real(ifft2(ifftshift(fxi1)))*Nx*Ny;
        fxi2 = fxi1(:,Nx/4+1:3*Nx/4);
        samp2 = real(ifft2(ifftshift(fxi2)))*Nx*Ny/2;
        varargout{1} = mu + samp1;
        varargout{2} = mu + samp2;
    case 3 % increments of increments along y
        fxi1 = sqrtcov1.*frr;
        samp1 = real(ifft2(ifftshift(fxi1)))*Nx*Ny;
        fxi2 = fxi1(Ny/4+1:3*Ny/4,:);
        samp2 = real(ifft2(ifftshift(fxi2)))*Nx*Ny/2;
        varargout{1} = mu + samp1;
        varargout{2} = mu + samp2;
    case 4 % increments of increments
        fxi1 = sqrtcov1.*frr;
        samp1 = real(ifft2(ifftshift(fxi1)))*Nx*Ny;
        fxi2 = fxi1(:,Nx/4+1:3*Nx/4);
        samp2 = real(ifft2(ifftshift(fxi2)))*Nx*Ny/2;
        fxi3 = fxi1(Ny/4+1:3*Ny/4,:);
        samp3 = real(ifft2(ifftshift(fxi3)))*Nx*Ny/2;
        fxi4 = fxi1(Ny/4+1:3*Ny/4,Nx/4+1:3*Nx/4);
        samp4 = real(ifft2(ifftshift(fxi4)))*Nx*Ny/4;
        varargout{1} = mu + samp1;
        varargout{2} = mu + samp2;
        varargout{3} = mu + samp3;
        varargout{4} = mu + samp4;
    case 5 % increments
        fxi1 = sqrtcov1.*frr;
        samp1 = real(ifft2(ifftshift(fxi1)))*Nx*Ny;
        fxi2 = fxi1(Ny/4+1:3*Ny/4,Nx/4+1:3*Nx/4);
        samp2 = real(ifft2(ifftshift(fxi2)))*Nx*Ny/4;
        varargout{1} = mu + samp1;
        varargout{2} = mu + samp2;
    otherwise
        fprintf(1,'Wrong Prior Type')
end

end

