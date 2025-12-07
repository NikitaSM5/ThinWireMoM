classdef AntennaGraph
    properties
        Vertices double = double.empty(0,3)   % координаты вершин графа
        Edges int32 = int32.empty(0,2)        % активные ребра
        AllowedEdges int32 = int32.empty(0,2) % все возможные ребра
        FeedEdge int32 = int32.empty(0,2)     % питающее ребро

        % Новые поля:
        AddableEdges int32 = int32.empty(0,2)    % ребра, которые можно добавить
        RemovableEdges int32 = int32.empty(0,2)  % ребра, которые можно удалить

        % --- Метаданные ---
        Meta struct = struct('ShapeType','','PlacementType','', ...
                             'GridType','', 'EdgeLength',NaN, 'Description','')

        ObjectiveValue double = NaN
    end

    methods
        function obj = AntennaGraph(vertices, edges, allowedEdges, meta, feedEdge)
            % Конструктор класса AntennaGraph

            if nargin == 0
                return;
            end

            if size(vertices,2) ~= 3
                error('Vertices должны быть размером [N x 3].');
            end

            if ~isempty(edges)
                edges = sort(edges, 2);
                edges = unique(edges, 'rows');
                edges(edges(:,1) == edges(:,2), :) = [];
            end

            if ~isempty(allowedEdges)
                allowedEdges = sort(allowedEdges, 2);
                allowedEdges = unique(allowedEdges, 'rows');
                allowedEdges(allowedEdges(:,1) == allowedEdges(:,2), :) = [];
            end

            if nargin < 4 || isempty(meta)
                meta = struct('ShapeType','','PlacementType','', ...
                              'GridType','', 'EdgeLength',NaN, 'Description','');
            end

            if nargin < 5 || isempty(feedEdge)
                feedEdge = int32.empty(0,2);
            else
                feedEdge = sort(feedEdge);
            end

            obj.Vertices     = vertices;
            obj.Edges        = edges;
            obj.AllowedEdges = allowedEdges;
            obj.Meta         = meta;
            obj.FeedEdge     = feedEdge;

            obj.ObjectiveValue = NaN;

            if ~isempty(edges) && ~isempty(allowedEdges)
                [tf, ~] = ismember(edges, allowedEdges, 'rows');
                if any(~tf)
                    warning('Некоторые активные рёбра не входят в список допустимых.');
                end
            end

            % --- первичный расчёт списков добавляемых/удаляемых рёбер ---
            obj.AddableEdges   = obj.getAddableEdges();
            obj.RemovableEdges = obj.getRemovableEdges();

            % --- опциональная проверка связности (как было) ---
            try
                if ~obj.isConnectedFromFeed() && ~isempty(obj.Edges)
                    warning('Созданный граф несвязен по активным рёбрам.');
                end
            catch
            end
        end


        function A = adjacencyMatrix(obj)
            N = size(obj.Vertices, 1);
            if isempty(obj.Edges)
                A = sparse(N, N);
                return;
            end
            i = obj.Edges(:,1);
            j = obj.Edges(:,2);
            v = ones(size(i));
            A = sparse([i; j], [j; i], [v; v], N, N);
            A(1:N+1:end) = 0;
        end


        function tf = isConnectedFromFeed(obj)
            if isempty(obj.Edges)
                tf = true;
                return;
            end

            if isempty(obj.FeedEdge)
                warning('FeedEdge не задан. Проверка связности невозможна.');
                tf = false;
                return;
            end

            activeVerts = unique(obj.Edges(:));
            A = obj.adjacencyMatrix();
            G = graph(A);
            bins = conncomp(G);

            feedNodes = obj.FeedEdge(:);
            feedNodes = feedNodes(feedNodes <= numel(bins));
            feedComp = bins(feedNodes(1));

            tf = all(bins(activeVerts) == feedComp);
        end


        function idx = findEdgeIndex(obj, edge)
            if isempty(obj.Edges)
                idx = [];
                return;
            end
            if numel(edge) ~= 2
                error('Ребро должно задаваться как [i j].');
            end
            edge = sort(edge);
            E = sort(obj.Edges, 2);
            [tf, loc] = ismember(edge, E, 'rows');
            if tf
                idx = loc;
            else
                idx = [];
            end
        end


        function tf = hasEdge(obj, edge)
            idx = obj.findEdgeIndex(edge);
            tf = ~isempty(idx);
        end


        function obj2 = addEdge(obj, edge)
            % Добавляет ребро и обновляет Addable/RemovableEdges

            obj2 = obj;

            if numel(edge) ~= 2
                error('Ребро должно задаваться как [i j].');
            end
            edge = sort(edge);

            if edge(1) == edge(2)
                warning('Петли не допускаются: [%d %d]', edge);
                return;
            end

            if ~isempty(obj.AllowedEdges)
                [isAllowed,~] = ismember(edge, sort(obj.AllowedEdges,2), 'rows');
                if ~isAllowed
                    warning('Ребро [%d %d] не входит в список допустимых.', edge);
                    return;
                end
            end

            if obj.hasEdge(edge)
                warning('Ребро [%d %d] уже добавлено.', edge);
                return;
            end

            obj2.Edges = [obj.Edges; edge];

            % Обновляем множества возможных рёбер
            obj2.AddableEdges   = obj2.getAddableEdges();
            obj2.RemovableEdges = obj2.getRemovableEdges();
        end


        function obj2 = removeEdge(obj, edge)
            % Удаляет ребро и обновляет Addable/RemovableEdges

            obj2 = obj;

            if numel(edge) ~= 2
                error('Ребро должно задаваться как [i j].');
            end
            edge = sort(edge);

            idx = obj.findEdgeIndex(edge);
            if isempty(idx)
                warning('Ребро [%d %d] не является активным. Отменено.', edge);
                return;
            end

            obj2.Edges(idx,:) = [];

            if ~isempty(obj.FeedEdge)
                if ~obj2.isConnectedFromFeed()
                    obj2 = obj;
                    warning('Удаление ребра [%d %d] нарушает связь с питанием. Отменено.', edge(1), edge(2));
                    return;
                end
            end

            % Обновляем множества возможных рёбер
            obj2.AddableEdges   = obj2.getAddableEdges();
            obj2.RemovableEdges = obj2.getRemovableEdges();
        end


        function addableEdges = getAddableEdges(obj)
            % Возвращает список рёбер, которые можно добавить в граф.
            %
            % Условия:
            %   - ребро должно быть в списке AllowedEdges
            %   - ребра ещё нет в графе (неактивное)
            %   - не петля (i ~= j)

            % --- Если список допустимых рёбер не задан ---
            if isempty(obj.AllowedEdges)
                addableEdges = int32.empty(0,2);
                return;
            end

            % --- Базовый набор кандидатов ---
            allowed = sort(int32(obj.AllowedEdges), 2);

            % Если активных рёбер нет — можно добавить все допустимые
            if isempty(obj.Edges)
                addableEdges = allowed;
                % на всякий случай убираем петли (их быть не должно)
                addableEdges(addableEdges(:,1) == addableEdges(:,2), :) = [];
                return;
            end

            % --- Есть активные рёбра: исключаем уже существующие ---
            existing = sort(int32(obj.Edges), 2);
            [~, idxMissing] = setdiff(allowed, existing, 'rows');

            addableEdges = allowed(idxMissing, :);
            addableEdges(addableEdges(:,1) == addableEdges(:,2), :) = [];
        end



        function removableEdges = getRemovableEdges(obj)
            % Рёбра, удаление которых не ломает связность относительно FeedEdge.
            % Если FeedEdge не задан — считаем, что ограничений нет, и
            % допускаем удаление любого ребра.

            if isempty(obj.Edges)
                removableEdges = int32.empty(0,2);
                return;
            end

            if isempty(obj.FeedEdge)
                removableEdges = obj.Edges;
                return;
            end

            removableEdges = int32.empty(0,2);

            for k = 1:size(obj.Edges, 1)
                edge = obj.Edges(k,:);

                tempGraph      = obj;
                tempGraph.Edges(k,:) = [];

                if tempGraph.isConnectedFromFeed()
                    removableEdges = [removableEdges; edge];
                end
            end
        end


        function obj2 = compress(obj)
            % Сжимает граф: оставляет только используемые вершины,
            % перенумеровывает рёбра и обновляет списки Addable/RemovableEdges.

            obj2 = obj;

            if isempty(obj.Edges)
                obj2.Vertices = obj.Vertices;
                return;
            end

            usedVerts = unique(obj.Edges(:));
            nUsed     = numel(usedVerts);

            map = zeros(size(obj.Vertices,1),1);
            map(usedVerts) = 1:nUsed;

            obj2.Vertices = obj.Vertices(usedVerts, :);
            obj2.Edges    = sort(map(obj.Edges), 2);

            if ~isempty(obj.AllowedEdges)
                validMask       = all(ismember(obj.AllowedEdges, usedVerts), 2);
                allowedFiltered = obj.AllowedEdges(validMask, :);
                obj2.AllowedEdges = sort(map(allowedFiltered), 2);
            end

            if ~isempty(obj.FeedEdge)
                if all(ismember(obj.FeedEdge, usedVerts))
                    obj2.FeedEdge = sort(map(obj.FeedEdge));
                else
                    warning('FeedEdge ссылался на удалённые вершины — сброшен.');
                    obj2.FeedEdge = int32.empty(0,2);
                end
            end

            obj2.Meta           = obj.Meta;
            obj2.ObjectiveValue = obj.ObjectiveValue;

            % Пересчёт множеств возможных рёбер
            obj2.AddableEdges   = obj2.getAddableEdges();
            obj2.RemovableEdges = obj2.getRemovableEdges();
        end


        function edgeCoords = getEdgeCoordinates(obj)
            M = size(obj.Edges, 1);
            if M == 0
                edgeCoords = zeros(0, 2, 3);
                return;
            end

            edgeCoords = zeros(M, 2, 3);
            V = obj.Vertices;

            for e = 1:M
                i = obj.Edges(e, 1);
                j = obj.Edges(e, 2);

                edgeCoords(e,1,1:3) = V(i,1:3);
                edgeCoords(e,2,1:3) = V(j,1:3);
            end
        end


        function [p1, p2] = getFeedEdgeCoordinates(obj)
            if isempty(obj.FeedEdge) || numel(obj.FeedEdge) ~= 2
                p1 = [];
                p2 = [];
                warning('FeedEdge не задан или имеет некорректный формат.');
                return;
            end

            i = obj.FeedEdge(1);
            j = obj.FeedEdge(2);

            nVerts = size(obj.Vertices, 1);
            if i < 1 || i > nVerts || j < 1 || j > nVerts
                p1 = [];
                p2 = [];
                warning('FeedEdge ссылается на несуществующие вершины.');
                return;
            end

            p1 = obj.Vertices(i, :);
            p2 = obj.Vertices(j, :);
        end
    end


    %%% Nikita`s Methods

    

    methods (Static)
        function obj = fromTemplate(templateParams)
            if nargin < 1 || isempty(templateParams)
                error('Не заданы параметры шаблона.');
            end

            [vertices, allowedEdges, meta] = GraphTemplates.generate(templateParams);

            if isfield(templateParams,'FeedEdge')
                feedEdge = templateParams.FeedEdge;
            else
                feedEdge = [];
            end

            obj = AntennaGraph(vertices, [], allowedEdges, meta, feedEdge);
        end

        function obj = emptyFromVertices(vertices, allowedEdges, meta)
            obj = AntennaGraph(vertices, [], allowedEdges, meta, []);
        end

        function obj = fullFromVertices(vertices, allowedEdges, meta)
            obj = AntennaGraph(vertices, allowedEdges, allowedEdges, meta, []);
        end
    end
end
