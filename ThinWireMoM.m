classdef ThinWireMoM < handle
    properties
        AntGraph  AntennaGraph = AntennaGraph()

        freq       double = NaN
        omega      double = NaN

        eps0       double = 8.854187817e-12
        mu0        double = 4*pi*1e-7
        eps_r      double = 1.0
        mu_r       double = 1.0

        eps        double = NaN
        mu         double = NaN

        k          double = NaN
        eta        double = NaN

        Segments   struct = struct([])
        Basis      struct = struct('Ndof',0,'vertexIds',[],'dofOfVertex',[],'dofs',[])

        WireRadius double = 1e-4

        Z   double = double.empty(0,0)
        V   double = double.empty(0,1)
        I   double = double.empty(0,1)
        NumIntPoints (1,1) double = 4

    end

    methods
        function obj = ThinWireMoM(antGraph, freq, eps_r, mu_r, wireRadius)
            obj.WireRadius = wireRadius;
            obj.AntGraph = antGraph;
            obj = obj.setMedium(eps_r, mu_r);
            obj = obj.setFrequency(freq);
        end

        function obj = setGraph(obj, antGraph)
            obj.AntGraph = antGraph;
        end

        function obj = setFrequency(obj, freqHz)
            obj.freq  = freqHz;
            obj.omega = 2*pi*freqHz;
            obj = obj.updateEMParams();
        end

        function obj = setMedium(obj, eps_r, mu_r)
            obj.eps_r = eps_r;
            obj.mu_r  = mu_r;
            obj.eps = obj.eps0 * eps_r;
            obj.mu  = obj.mu0  * mu_r;
            obj = obj.updateEMParams();
        end

        function lambda = getWavelength(obj)
                lambda = 2*pi / obj.k;
        end

        function k = getWaveNumber(obj)
            k = obj.k;
        end

        function eta = getWaveImpedance(obj)
            eta = obj.eta;
        end

        function obj = prepareGeometryAndBasis(obj)
            obj = obj.buildSegmentsFromGraph();
            obj = obj.buildVertexHatBasis();
        end

        function obj = buildSegmentsFromGraph(obj)
            E = obj.AntGraph.Edges;
            V = obj.AntGraph.Vertices;

            M = size(E,1);
            Seg(M,1) = struct('id',0,'v1',int32(0),'v2',int32(0), ...
                'r1',[],'r2',[],'center',[],'length',0.0,'t',[]);

            for s = 1:M
                v1 = E(s,1);
                v2 = E(s,2);

                r1 = V(v1, :);
                r2 = V(v2, :);

                d  = r2 - r1;
                L  = norm(d);

                t = d / L;

                Seg(s).id     = s;
                Seg(s).v1     = int32(v1);
                Seg(s).v2     = int32(v2);
                Seg(s).r1     = r1;
                Seg(s).r2     = r2;
                Seg(s).center = 0.5*(r1 + r2);
                Seg(s).length = L;
                Seg(s).t      = t;
            end

            obj.Segments = Seg;
        end

        function obj = buildVertexHatBasis(obj)
            Seg = obj.Segments;
            E   = obj.AntGraph.Edges;

            activeVerts = unique(E(:));
            Ndof        = numel(activeVerts);

            nVerts       = size(obj.AntGraph.Vertices, 1);
            dofOfVertex  = zeros(nVerts,1,'int32');

            for k = 1:Ndof
                v = activeVerts(k);
                dofOfVertex(v) = int32(k);
            end

            M = numel(Seg);
            vertexSegList = cell(nVerts,1);
            for v = 1:nVerts
                vertexSegList{v} = int32([]);
            end

            for s = 1:M
                v1 = Seg(s).v1;
                v2 = Seg(s).v2;
                vertexSegList{v1}(end+1,1) = int32(s);
                vertexSegList{v2}(end+1,1) = int32(s);
            end


            dofs(1,Ndof) = struct('vertexId',int32(0), 'incidentSegs', []);

            for k = 1:Ndof
                v = activeVerts(k);
                segsForVertex = vertexSegList{v};

                incSegs = struct('segId',{}, 'sign',{}, 'length',{});

                for iseg = 1:numel(segsForVertex)
                    sId = segsForVertex(iseg);
                    s   = Seg(sId);

                    if s.v1 == v
                        signLocal = +1;
                    elseif s.v2 == v
                        signLocal = -1;
                    end

                    incSegs(end+1).segId  = int32(sId);
                    incSegs(end).sign     = signLocal;
                    incSegs(end).length   = s.length;
                end

                dofs(k).vertexId     = int32(v);
                dofs(k).incidentSegs = incSegs;
            end

            Basis.Ndof        = Ndof;
            Basis.vertexIds   = int32(activeVerts(:));
            Basis.dofOfVertex = dofOfVertex;
            Basis.dofs        = dofs;

            obj.Basis = Basis;
        end

        function [Lambda, dLambda_ds] = evalHatOnSegment(obj, segId, s_local)
            seg = obj.Segments(segId);
            L   = seg.length;

            Lambda      = s_local / L;
            dLambda_ds  = 1 / L;
        end

        function t_local = getLocalTangent(obj, dofIdx, segId)
            s   = obj.Segments(segId);
            v   = obj.Basis.dofs(dofIdx).vertexId;

            % if s.v1 == v
            %     signLocal = +1;
            % elseif s.v2 == v
            %     signLocal = -1;
            % end

            % t_local = signLocal * s.t;
            t_local =  s.t;
        end


        function obj = assembleZ(obj, numIntPoints)

            obj.NumIntPoints = numIntPoints;
            Ndof = obj.Basis.Ndof;

            Zmat = complex(zeros(Ndof, Ndof));

            for m = 1:Ndof
                dof_m = obj.Basis.dofs(m).incidentSegs;
                for n = 1:Ndof
                    dof_n = obj.Basis.dofs(n).incidentSegs;

                    Zmn = 0;

                    for im = 1:numel(dof_m)
                        for jn = 1:numel(dof_n)
                            Zmn = Zmn + obj.segmentPairContribution(m, n, ...
                                dof_m(im), dof_n(jn));
                        end
                    end

                    Zmat(m,n) = Zmn;
                end
            end

            obj.Z = Zmat;
        end

        function obj = assembleVDeltaGap(obj, Vgap)
            Ndof = obj.Basis.Ndof;

            feedEdge = obj.AntGraph.FeedEdge;
            feedEdge = sort(double(feedEdge(:))).';

            segId = obj.findFeedSegmentId(feedEdge);

            dofOfVertex = obj.Basis.dofOfVertex;
            v1 = feedEdge(1);
            v2 = feedEdge(2);

            dof1 = dofOfVertex(v1);
            dof2 = dofOfVertex(v2);

            Vvec = complex(zeros(Ndof,1));

            incSegs1 = obj.Basis.dofs(dof1).incidentSegs;
            for k = 1:numel(incSegs1)
                if incSegs1(k).segId == segId
                    signLocal = double(incSegs1(k).sign);
                    Vvec(dof1) = Vvec(dof1) + signLocal * (Vgap/2);
                    break;
                end
            end

            incSegs2 = obj.Basis.dofs(dof2).incidentSegs;
            for k = 1:numel(incSegs2)
                if incSegs2(k).segId == segId
                    signLocal = double(incSegs2(k).sign);
                    Vvec(dof2) = Vvec(dof2) + signLocal * (Vgap/2);
                    break;
                end
            end

            obj.V = Vvec;
        end

        function segId = findFeedSegmentId(obj, feedEdge)
            feedEdge = sort(double(feedEdge(:))).';

            Seg = obj.Segments;
            segId = [];

            for s = 1:numel(Seg)
                v1v2 = sort(double([Seg(s).v1, Seg(s).v2]));
                if v1v2(1) == feedEdge(1) && v1v2(2) == feedEdge(2)
                    segId = s;
                    break;
                end
            end
        end


    end

    methods (Access = private)
        function obj = updateEMParams(obj)
            obj.k   = obj.omega * sqrt(obj.mu * obj.eps);
            obj.eta = sqrt(obj.mu / obj.eps);
        end

        function Zmn = segmentPairContribution(obj, m, n, segInfo_m, segInfo_n)

            seg_m = obj.Segments(segInfo_m.segId);
            seg_n = obj.Segments(segInfo_n.segId);

            Lm = segInfo_m.length;
            Ln = segInfo_n.length;

            sign_m = double(segInfo_m.sign);
            sign_n = double(segInfo_n.sign);

            isSameDof    = (m == n);
            isSameSegment = (segInfo_m.segId == segInfo_n.segId);

            if isSameDof && isSameSegment
                % Zmn = obj.selfTermSameSegment(seg_m, segInfo_m);
                % return;
            end

            Nq = obj.NumIntPoints;

            ds_m = Lm / Nq;
            ds_n = Ln / Nq;

            t_m = obj.getLocalTangent(m, segInfo_m.segId);
            t_n = obj.getLocalTangent(n, segInfo_n.segId);
            tdot = dot(t_m, t_n);

            sumJ = 0;
            sumD = 0;

            for im = 1:Nq
                s_m = (im - 0.5) * ds_m;

                r_m = obj.evalPointOnSegmentLocal(seg_m, s_m, sign_m);

                [Lambda_m, dLambda_m_local] = obj.evalHatOnSegment(segInfo_m.segId, s_m);

                dLambda_m = sign_m * dLambda_m_local;

                for jn = 1:Nq
                    s_n = (jn - 0.5) * ds_n;
                    r_n = obj.evalPointOnSegmentLocal(seg_n, s_n, sign_n);

                    Rvec = r_m - r_n;
                    R    = norm(Rvec);

                    [Lambda_n, dLambda_n_local] = obj.evalHatOnSegment(segInfo_n.segId, s_n);
                    dLambda_n = sign_n * dLambda_n_local;

                    G = obj.greenFunction(R);

                    sumJ = sumJ + Lambda_m * Lambda_n * tdot * G;
                    sumD = sumD + dLambda_m * dLambda_n * G;
                end
            end
            Zmn =  (sumJ - 1/obj.k^2 * sumD);
            
        end

        function G = greenFunction(obj, R)
            a    = obj.WireRadius;
            Reff = sqrt(R.^2 + a.^2);
            G    = exp(-1i * obj.k * Reff) ./ (4*pi*Reff);
        end



        function r = evalPointOnSegmentLocal(~, seg, s_local, signLocal)
            if signLocal > 0
                r = seg.r1 + seg.t * s_local;
            else
                r = seg.r2 - seg.t * s_local;
            end
        end


        function S1 = S1_triangle(obj, x, Dl)
            a = obj.WireRadius;
            k = obj.k;

            sqrt1 = sqrt(a^2 + (x - Dl).^2);
            sqrt2 = sqrt(a^2 + x.^2);

            num = x + sqrt2;
            den = x - Dl + sqrt1;

            S1 = (sqrt1 - sqrt2)/Dl ...
                + (x./Dl).*log(num./den) ...
                - 1i * k * Dl / 2;
        end

        function S2 = S2_triangle(obj, x, Dl, sameSign)
            a = obj.WireRadius;
            k = obj.k;

            sqrt1 = sqrt(a^2 + (x - Dl).^2);
            sqrt2 = sqrt(a^2 + x.^2);

            num = x + sqrt2;
            den = x - Dl + sqrt1;
            logTerm = log(num./den);

            sgn = 1;
            if ~sameSign
                sgn = -1;
            end

            S2 = sgn * (logTerm - 1i * k * Dl) / (Dl^2);
        end

        function Zself = selfTermSameSegment(obj, seg, segInfo)

            k   = obj.k;
            L   = segInfo.length;
            sgnVert = double(segInfo.sign);

            Nq = obj.NumIntPoints;

            dx = L / Nq;

            acc = 0;

            for p = 1:Nq

                x = (p - 0.5) * dx;

                % fm = obj.evalHatOnSegment(segInfo.segId, x);

                S1 = obj.S1_triangle(x, L);
                S2 = obj.S2_triangle(x, L, true);

                acc = acc + ( fm * S1 - (1/k^2) * S2 );
            end

            Zself =  (1/(4*pi)) * acc * dx;
        end

    end
end
