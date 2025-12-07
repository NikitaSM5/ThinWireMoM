classdef GraphVisualizer
    %GRAPHVISUALIZER Визуализация объектов AntennaGraph
    %
    %   h = GraphVisualizer.plotGraph(G, Name,Value,...)
    %
    %   Имя–значение:
    %       'Title'          - заголовок фигуры (string), по умолчанию авто
    %       'ShowNodeLabels' - логическое, подписывать ли номера вершин (false)
    %
    %   Возвращает хэндл фигуры h.
    %
    %   GraphVisualizer.saveFigure(h, 'basename')
    %       сохраняет фигуру в basename.png и basename.eps

    methods (Static)
        function hFig = plotGraph(G, varargin)
            % Проверка входа
            if ~isa(G, 'AntennaGraph')
                error('plotGraph ожидает объект AntennaGraph.');
            end

            % ---- Парсинг опций ----
            p = inputParser;
            addParameter(p, 'Title', '', @(x)ischar(x) || isstring(x));
            addParameter(p, 'ShowNodeLabels', false, @(x)islogical(x) || isnumeric(x));
            parse(p, varargin{:});
            opt = p.Results;

            % ---- Геометрия ----
            V = G.Vertices;
            if isempty(V)
                error('У графа нет вершин.');
            end

            % Решаем, 2D или 3D
            z = V(:,3);
            if max(abs(z - z(1))) < 1e-9
                is3D = false;   % все z примерно одинаковы -> плоский граф
            else
                is3D = true;
            end

            % ---- Предрасчёт множеств рёбер ----
            activeEdges    = G.Edges;
            feedEdge       = G.FeedEdge;
            addableEdges   = G.getAddableEdges();
            removableEdges = G.getRemovableEdges();

            % ---- Фигура и оси ----
            hFig = figure('Color','w','Name','Antenna graph');
            hAx  = axes('Parent',hFig);
            hold(hAx, 'on'); grid(hAx,'on'); box(hAx,'on');

            if is3D
                view(hAx, 3);
            else
                view(hAx, 2);
            end

            axis(hAx, 'equal');

            % -----------------------------------------
            %   Рисуем вершины
            % -----------------------------------------
            if is3D
                hNodes = plot3(hAx, V(:,1), V(:,2), V(:,3), ...
                    'o', 'MarkerFaceColor',[0.8 0.8 0.8], ...
                    'MarkerEdgeColor',[0.3 0.3 0.3], ...
                    'MarkerSize',5);
            else
                hNodes = plot(hAx, V(:,1), V(:,2), ...
                    'o', 'MarkerFaceColor',[0.8 0.8 0.8], ...
                    'MarkerEdgeColor',[0.3 0.3 0.3], ...
                    'MarkerSize',5);
            end

            % Подписи вершин (по желанию)
            if opt.ShowNodeLabels
                for v = 1:size(V,1)
                    if is3D
                        text(hAx, V(v,1), V(v,2), V(v,3)+0.01, sprintf('%d',v), ...
                            'FontSize',8,'HorizontalAlignment','center');
                    else
                        text(hAx, V(v,1), V(v,2), sprintf('%d',v), ...
                            'FontSize',8,'HorizontalAlignment','center', ...
                            'VerticalAlignment','bottom');
                    end
                end
            end

            % -----------------------------------------
            %   Вспомогательная функция рисования рёбер
            % -----------------------------------------
            function hLines = drawEdges(E, color, style, width, visible)
                if isempty(E)
                    hLines = gobjects(0,1);
                    return;
                end
                nE = size(E,1);
                hLines = gobjects(nE,1);
                for ee = 1:nE
                    i = E(ee,1);
                    j = E(ee,2);
                    if is3D
                        hLines(ee) = plot3(hAx, ...
                            [V(i,1) V(j,1)], [V(i,2) V(j,2)], [V(i,3) V(j,3)], ...
                            'Color', color, 'LineStyle',style, 'LineWidth',width, ...
                            'Visible', visible);
                    else
                        hLines(ee) = plot(hAx, ...
                            [V(i,1) V(j,1)], [V(i,2) V(j,2)], ...
                            'Color', color, 'LineStyle',style, 'LineWidth',width, ...
                            'Visible', visible);
                    end
                end
            end

            % -----------------------------------------
            %   Активные рёбра (синие)
            % -----------------------------------------
            hActive = drawEdges(activeEdges, [0 0 0.8], '-', 1.5, 'on');

            % -----------------------------------------
            %   Питающее ребро (фиолетовое, толстое)
            % -----------------------------------------
            if ~isempty(feedEdge)
                hFeed = drawEdges(feedEdge(:).', [0.7 0 0.7], '-', 3, 'on');
                % Выделим вершины питания другим цветом
                vFeed = unique(feedEdge(:));
                if ~isempty(vFeed)
                    if is3D
                        plot3(hAx, V(vFeed,1), V(vFeed,2), V(vFeed,3), ...
                            'o', 'MarkerSize',7, ...
                            'MarkerFaceColor',[1 0.6 0], ...
                            'MarkerEdgeColor',[0.7 0.3 0]);
                    else
                        plot(hAx, V(vFeed,1), V(vFeed,2), ...
                            'o', 'MarkerSize',7, ...
                            'MarkerFaceColor',[1 0.6 0], ...
                            'MarkerEdgeColor',[0.7 0.3 0]);
                    end
                end
            else
                hFeed = gobjects(0,1);
            end

            % -----------------------------------------
            %   Добавляемые рёбра (зелёные, по умолчанию скрыты)
            % -----------------------------------------
            hAddable = drawEdges(addableEdges, [0 0.6 0], '-', 1, 'off');

            % -----------------------------------------
            %   Удаляемые рёбра (красные, по умолчанию скрыты)
            % -----------------------------------------
            hRemovable = drawEdges(removableEdges, [0.9 0 0], '-', 1.5, 'off');

            % Подписи осей
            xlabel(hAx, 'X'); ylabel(hAx, 'Y');
            if is3D, zlabel(hAx,'Z'); end

            % Заголовок
            if ~isempty(opt.Title)
                title(hAx, opt.Title, 'FontWeight','bold');
            else
                if is3D
                    tdim = '3D';
                else
                    tdim = '2D';
                end
                title(hAx, sprintf('AntennaGraph (%s, %s grid)', ...
                    tdim, G.Meta.GridType), 'FontWeight','bold');
            end

            % -----------------------------------------
            %   Чекбоксы управления слоями
            % -----------------------------------------
            % Позиции в нормированных координатах
            cbWidth  = 0.18;
            cbHeight = 0.04;

            % Добавляемые
            hCBadd = uicontrol('Style','checkbox',...
                'String','Show addable edges (green)', ...
                'Units','normalized', ...
                'Position',[0.01 0.01 cbWidth cbHeight], ...
                'Value',0, ...
                'Callback',@toggleAddable);

            if isempty(hAddable)
                set(hCBadd,'Enable','off');
            end

            % Удаляемые
            hCBrem = uicontrol('Style','checkbox',...
                'String','Show removable edges (red)', ...
                'Units','normalized', ...
                'Position',[0.25 0.01 cbWidth cbHeight], ...
                'Value',0, ...
                'Callback',@toggleRemovable);

            if isempty(hRemovable)
                set(hCBrem,'Enable','off');
            end

            % Вспомогательные callback-и
            function toggleAddable(src,~)
                if get(src,'Value')
                    set(hAddable,'Visible','on');
                else
                    set(hAddable,'Visible','off');
                end
            end

            function toggleRemovable(src,~)
                if get(src,'Value')
                    set(hRemovable,'Visible','on');
                else
                    set(hRemovable,'Visible','off');
                end
            end
        end

        function saveFigure(hFig, baseFilename)
            %SAVEFIGURE Сохранить фигуру в PNG и EPS
            %
            %   GraphVisualizer.saveFigure(hFig, 'my_graph')
            %   создаст файлы my_graph.png и my_graph.eps

            if nargin < 1 || isempty(hFig)
                hFig = gcf;
            end
            if nargin < 2 || isempty(baseFilename)
                baseFilename = 'graph';
            end

            if ~ishandle(hFig) || ~strcmp(get(hFig,'Type'),'figure')
                error('Первый аргумент должен быть хэндлом фигуры.');
            end

            pngFile = [baseFilename, '.png'];
            epsFile = [baseFilename, '.eps'];

            % Можно настроить разрешение при желании
            print(hFig, pngFile, '-dpng', '-r300');
            print(hFig, epsFile, '-depsc');

            fprintf('Figure saved to %s and %s\n', pngFile, epsFile);
        end
    end
end
