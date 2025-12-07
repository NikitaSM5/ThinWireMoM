clc; clear;
close all;

eps_r = 1.0;
mu_r  = 1.0;

c = 3e8;  
f = 300e6;          
lambda = c / f;     
a = lambda/1e9;
L_factor = 3;       
L        = L_factor*lambda*1;         

Nseg = 101;
fprintf('f = %.1f МГц, λ = %.4f м, L = %.4f м (%.2f λ), Nseg = %d\n', ...
    f/1e6, lambda, L, L/lambda, Nseg);

Nz = Nseg + 1;
z = linspace(-L/2, L/2, Nz).';
Vertices = [zeros(Nz,2), z]; 

Edges = int32([ (1:Nseg).', (2:Nseg+1).' ]);

allowedEdges = Edges;

v_feed1 = floor(Nseg/2) + 1;
v_feed2 = floor(Nseg/2) + 2;
FeedEdge = int32([v_feed1, v_feed2]);

meta = struct('ShapeType','dipole', ...
              'PlacementType','', ...
              'GridType','', ...
              'EdgeLength',L/Nseg, ...
              'Description','Straight center-fed dipole for ThinWireMoM test');

G = AntennaGraph(Vertices, Edges, allowedEdges, meta, FeedEdge);
GraphVisualizer.plotGraph(G, 'Title','dipole', 'ShowNodeLabels',false);

mom = ThinWireMoM(G, f, eps_r, mu_r, a);

mom = mom.prepareGeometryAndBasis();

NumIntPoints = 8;
mom = mom.assembleZ(NumIntPoints);

Vgap = 1.0;
mom = mom.assembleVDeltaGap(Vgap);

mom.I = mom.Z \ mom.V;

segId = mom.findFeedSegmentId(FeedEdge);

dofOfVertex = mom.Basis.dofOfVertex;
dof1 = dofOfVertex(v_feed1);
dof2 = dofOfVertex(v_feed2);


M = numel(mom.Segments);
Iseg = zeros(M,1);

for s = 1:M
    seg = mom.Segments(s);

    v1 = double(seg.v1);
    v2 = double(seg.v2);
    d1 = dofOfVertex(v1);
    d2 = dofOfVertex(v2);

    Iseg(s) = 0.5*(mom.I(d1) - mom.I(d2));
end


Iseg(51) = (Iseg(50) + Iseg(52))/2;

Zin    = Vgap / Iseg(51);

fprintf('Z_in = %.3f %+.3fj Ом\n', real(Zin), imag(Zin));



figure;
plot(1:M, abs(Iseg/abs(max(Iseg))), '-o'); grid on; 
xlim([0, M]);
xlabel('Номер ребра');
ylabel('|I|');

hold on

z_seg = zeros(M,1);
for s = 1:M
    seg = mom.Segments(s);
    v1 = double(seg.v1);
    v2 = double(seg.v2);
    z_seg(s) = 0.5*(Vertices(v1,3) + Vertices(v2,3));
end



k0 = 2*pi/lambda;
I_theory = sin(k0*(L/2 - abs(z_seg)));

plot(1:M, real(I_theory), 'r-', 'LineWidth', 2);
legend('MoM','Теоретическое распределение','Location','best');
