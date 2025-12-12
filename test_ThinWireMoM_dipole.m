clc; clear;
close all;

eps_r = 1.0;
mu_r  = 1.0;
c = 299792458;  
f = 300e6;          
lambda = c / f;   
NumIntPoints = 8;
a = lambda/1e4;
L_factor = 0.47;       
L        = L_factor*lambda*1;     

Nseg = 31;
dL = L/Nseg/NumIntPoints;


fprintf('f = %.1f МГц, λ = %.4f м, L = %.4f м (%.2f λ), a = %.4f, dL = %.4f, dL/a = %.4f ,Nseg = %d\n', ...
    f/1e6, lambda, L, L/lambda, a , dL, dL/a , Nseg);

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


mom = mom.assembleZ(NumIntPoints);

Vgap = 1.0;
mom = mom.assembleVDeltaGap(Vgap);

mom.I = mom.Z \ mom.V;

nVerts = size(Vertices,1);
Ivert = zeros(nVerts,1);
for v = 1:nVerts-1
    d = mom.Basis.dofOfVertex(v);
    if d > 0
        Ivert(v) = mom.I(d);
    else
        Ivert(v) = 0;
    end
end

M = numel(mom.Segments);
Iseg = zeros(M,1);
for s = 1:M
    v1 = double(mom.Segments(s).v1);
    v2 = double(mom.Segments(s).v2);
    Iseg(s) = 0.5*(Ivert(v1) + Ivert(v2));
end

segId = mom.findFeedSegmentId(FeedEdge);
If = Iseg(segId);
Zin = Vgap / If;
fprintf('Zin = %.6f %+.6fj Ом\n', real(Zin), imag(Zin));
fprintf('Rin = %.6f Ом\n', abs(Zin));
figure;
plot(abs(Iseg), '-o'); grid on;
xlim([0, M]);
xlabel('Номер ребра');
ylabel('|I|');




