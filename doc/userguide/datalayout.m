
Ny = 8;
Nx = 12;

clf; 
hold on;

%axis([-3 Nx+2 -3 Ny+2])
%axis off
%axis equal
%axis ij

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real-valued data layout

% draw grey patch to indicate padding memory
light = [0.9 0.9 0.9];
patch([Nx-2, Nx, Nx, Nx-2], [0 0 Ny Ny], light)

tarrow(0.5, 0.5, 11.8, 0.5, .4, ' ');

% draw horizontal lines
for ny=0:Ny 
  plot([0 Nx], [ny ny], 'k')
end

% draw vertical lines, some dashed
for nx=0:Nx
  plot([nx nx], [0 Ny], 'k')
end

% add labels
for n=0:Nx-3;
  adjust = 0.35;
  %if n>=10 adjust =0.2; end;
  text(n+adjust, -0.5, int2str(n));
end
text(Nx+1, -0.5, 'nz');

% add labels
for ny=0:Ny-1;
  text(-.7, ny+0.5, int2str(ny));
end
text(-.8, Ny+0.5, 'nx');

text(1, Ny+0.5, 'real-valued data f(nx,nz)')
text(9.8, Ny+0.5, '(padding)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Complex-valued data layout

M=Ny+4; % y-offset for entire complex chart

% draw horizontal lines
for ny=0:Ny 
  plot([0 Nx], [M+ny M+ny], 'k')
end

% draw vertical lines, some dashed
for nx=0:2:Nx
  plot([nx nx], [M+0 M+Ny], 'k')
end
for nx=1:2:Nx-1
  plot([nx nx], [M+0 M+Ny], 'k:')
end

% add labels
for n=0:Nx/2-1;
  text(2*n+0.9, M-0.5, int2str(n));
  text(2*n+0.9, M-1.2, int2str(n));
end
text(Nx+1, M-1.2, 'mz');
text(Nx+1, M-0.5, 'kz');

% add labels
for ny=0:Ny-1;
  text(-1.5, M+ny+0.5, int2str(ny));
end
text(-1.6, M+Ny+0.6, 'mx');
for ny=0:Ny/2;
  text(-0.5, M+ny+0.5, int2str(ny));
end
for ny=Ny/2+1:Ny-1;
  text(-0.8, M+ny+0.5, int2str(-Ny+ny));
end
text(-0.6, M+Ny+0.6, 'kx');

text(2, M+Ny+0.5, 'complex-valued data f^{ \sim}(nx,nz)')


%axis([-3 Nx+2 -3 Ny+2])
axis off
axis equal
axis ij

print -deps2 datalayout.eps

