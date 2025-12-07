classdef GraphTemplates
    % Класс для генерации шаблонных сеток графов
    % Реализует цилиндр (surface), с вариантами rect и tri.
    %
    % Поддерживается авто-подбор:
    %   - NSides (число сегментов по окружности)
    %   - NLevels (число уровней по высоте)
    % по заданным Radius, Height и целевой EdgeLength.

    methods (Static)
        function [vertices, allowedEdges, meta] = generate(templateParams)
            % Генерация геометрии и связей для шаблонного графа.
            %
            % Обязательные поля:
            %   .ShapeType     = 'cylinder'
            %   .PlacementType = 'surface'
            %   .GridType      = 'rect' | 'tri'
            %   .Radius        - радиус цилиндра
            %   .Height        - высота цилиндра
            %
            % Варианты задания:
            %   1) Пользователь задаёт NSides и NLevels явно
            %   2) Пользователь задаёт EdgeLength, а NSides/NLevels
            %      автоматически подбираются по минимальной ошибке:
            %          - по окружности: длина хорды ~ EdgeLength
            %          - по высоте: шаг по Z ~ EdgeLength
            %
            % Выход:
            %   vertices     — [V x 3] координаты вершин
            %   allowedEdges — [E x 2] индексы допустимых рёбер
            %   meta         — структура описания

            if isfield(templateParams, 'ShapeType') && strcmpi(templateParams.ShapeType, 'rectangle')
                % =====================================================
                % === Прямоугольная область (Z=0)
                % =====================================================
                % Поля:
                %   .Width, .Height, .EdgeLength, .GridType = 'rect' | 'tri'
                % =====================================================

                if ~isfield(templateParams, 'GridType')
                    error('Для rectangle нужно указать GridType = "rect" или "tri".');
                end
                gridType = lower(templateParams.GridType);

                if ~isfield(templateParams, 'Width') || ~isfield(templateParams, 'Height')
                    error('Для rectangle необходимо задать Width и Height.');
                end
                W = templateParams.Width;
                H = templateParams.Height;

                if ~isfield(templateParams, 'EdgeLength') || isempty(templateParams.EdgeLength)
                    error('Для rectangle необходимо задать EdgeLength.');
                end
                L = templateParams.EdgeLength;

                % --- Количество узлов по каждой оси ---
                Nx = max(2, round(W / L) + 1);
                Ny = max(2, round(H / L) + 1);
                dx = W / (Nx - 1);
                dy = H / (Ny - 1);

                switch gridType
                    % =====================================================
                    % === Прямоугольная сетка (Манхэттен)
                    % =====================================================
                    case 'rect'
                        % Узлы прямоугольной решётки в "тетрадном" порядке:
                        % индексы идут по строкам: слева направо, снизу вверх.
                        vertices = zeros(Nx*Ny, 3);
                        k = 0;
                        for iy = 1:Ny
                            y = (iy-1)*dy;
                            for ix = 1:Nx
                                x = (ix-1)*dx;
                                k = k + 1;
                                vertices(k,:) = [x, y, 0];
                            end
                        end

                        allowedEdges = int32.empty(0,2);

                        % Связи только с ближайшими соседями:
                        % вправо и вверх (влево/вниз автоматически
                        % получаются, так как граф неориентированный)
                        for iy = 1:Ny
                            for ix = 1:Nx
                                idx = int32((iy-1)*Nx + ix);  % теперь это реально (строка, столбец)

                                % вправо: тот же ряд, следующий столбец
                                if ix < Nx
                                    allowedEdges(end+1,:) = [idx, idx+1];
                                end

                                % вверх: тот же столбец, следующий ряд
                                if iy < Ny
                                    allowedEdges(end+1,:) = [idx, idx+Nx];
                                end
                            end
                        end

                        allowedEdges = unique(sort(allowedEdges,2),'rows');


                        % =====================================================
                        % === Треугольная сетка (смещённые строки)
                        % =====================================================
                    case 'tri'
                        % Строим ряды, каждый второй смещён на dx/2
                        vertices = [];
                        for iy = 1:Ny
                            shift = mod(iy-1,2) * dx/2;      % 0, dx/2, 0, dx/2, ...
                            xRow  = 0:dx:W;
                            xRow  = xRow + shift;
                            xRow(xRow > W) = [];            % отсекаем вылезшие за границу
                            yRow  = (iy-1)*dy * ones(size(xRow));
                            vertices = [vertices; [xRow', yRow', zeros(numel(xRow),1)]];
                        end

                        allowedEdges = int32.empty(0,2);

                        % --- Горизонтальные связи внутри каждого ряда ---
                        for iy = 1:Ny
                            y  = (iy-1)*dy;
                            row = find(abs(vertices(:,2)-y) < 1e-9);
                            row = sort(row);
                            for k = 1:numel(row)-1
                                allowedEdges(end+1,:) = int32([row(k), row(k+1)]);
                            end
                        end

                        % --- Диагональные связи между соседними рядами ---
                        % Правило: соединяем вершины только с соседнего ряда,
                        % причём только те, для которых |x2 - x1| <= dx (иначе это
                        % прыжок через колонку).
                        for iy = 1:Ny-1
                            y1 = (iy-1)*dy;
                            y2 = iy*dy;
                            row1 = find(abs(vertices(:,2)-y1) < 1e-9);
                            row2 = find(abs(vertices(:,2)-y2) < 1e-9);
                            row1 = sort(row1);
                            row2 = sort(row2);

                            x2_all = vertices(row2,1);

                            for k = 1:numel(row1)
                                x1 = vertices(row1(k),1);
                                d  = abs(x2_all - x1);

                                % берём только соседей по x в пределах одного шага
                                neighMask = d <= dx*1.01;
                                neighIdx  = row2(neighMask);

                                % обычно 1 или 2 вершины, соединяем со всеми
                                for n = neighIdx'
                                    allowedEdges(end+1,:) = int32([row1(k), n]);
                                end
                            end
                        end

                        allowedEdges = unique(sort(allowedEdges,2),'rows');

                    otherwise
                        error('GridType должен быть "rect" или "tri".');
                end

                % --- Метаданные ---
                meta = struct( ...
                    'ShapeType', 'rectangle', ...
                    'GridType', gridType, ...
                    'EdgeLength', L, ...
                    'Description', sprintf('Rectangular surface %.1fx%.1f (%s grid)', W, H, gridType));

                return;
            end






            % --- Проверки типа фигуры ---
            if ~isfield(templateParams, 'ShapeType') || ~strcmpi(templateParams.ShapeType, 'cylinder')
                error('Сейчас поддерживается только ShapeType = "cylinder".');
            end
            if ~isfield(templateParams, 'PlacementType') || ~strcmpi(templateParams.PlacementType, 'surface')
                error('Сейчас поддерживается только PlacementType = "surface".');
            end

            if ~isfield(templateParams, 'GridType')
                error('Не указан GridType ("rect" или "tri").');
            end
            gridType = lower(templateParams.GridType);

            % --- Основные геометрические параметры ---
            R  = templateParams.Radius;
            H  = templateParams.Height;

            % --- Целевая длина ребра (опционально, но нужна для авто-подбора) ---
            L = NaN;
            if isfield(templateParams, 'EdgeLength') && ~isempty(templateParams.EdgeLength)
                L = templateParams.EdgeLength;
            end

            % --- NSides (по окружности) ---
            if isfield(templateParams, 'NSides') && ~isempty(templateParams.NSides)
                N = templateParams.NSides;
            else
                if isnan(L)
                    error('Для авто-подбора NSides требуется EdgeLength.');
                end
                N = GraphTemplates.estimateCylinderNSides(R, L);
            end

            % --- NLevels (по высоте) ---
            if isfield(templateParams, 'NLevels') && ~isempty(templateParams.NLevels)
                M = templateParams.NLevels;
            else
                if isnan(L)
                    error('Для авто-подбора NLevels требуется EdgeLength.');
                end
                M = GraphTemplates.estimateCylinderNLevels(H, L);
            end

            % --- Описание по умолчанию ---
            if ~isfield(templateParams, 'Description') || isempty(templateParams.Description)
                templateParams.Description = sprintf('Cylindrical surface (%s grid), N=%d, M=%d', gridType, N, M);
            end

            % --- Углы и уровни ---
            theta   = linspace(0, 2*pi, N+1); theta(end) = [];
            zLevels = linspace(0, H, M);

            % ===========================================================
            % Прямоугольная сетка
            % ===========================================================
            if strcmpi(gridType, 'rect')
                vertices = zeros(N*M, 3);

                for m = 1:M
                    th = theta;
                    x = R * cos(th);
                    y = R * sin(th);
                    z = zLevels(m) * ones(1, N);
                    vertices((m-1)*N+1:m*N,:) = [x' y' z'];
                end

                allowedEdges = [];

                % Горизонтальные связи по окружности
                for m = 1:M
                    base = (m-1)*N;
                    for k = 1:N
                        i = base + k;
                        j = base + mod(k,N) + 1;
                        allowedEdges = [allowedEdges; i j];
                    end
                end

                % Вертикальные связи
                for m = 1:M-1
                    base     = (m-1)*N;
                    nextBase = m*N;
                    for k = 1:N
                        i = base + k;
                        j = nextBase + k;
                        allowedEdges = [allowedEdges; i j];
                    end
                end

            % ===========================================================
            % Треугольная сетка
            % ===========================================================
            elseif strcmpi(gridType, 'tri')
                % Геометрия: каждый уровень поворачиваем на половину
                % углового шага относительно предыдущего.
                dtheta   = 2*pi / N;
                vertices = zeros(N*M, 3);

                for m = 1:M
                    th = theta + (m-1)*dtheta/2;
                    x  = R * cos(th);
                    y  = R * sin(th);
                    z  = zLevels(m) * ones(1,N);
                    vertices((m-1)*N+1:m*N,:) = [x' y' z'];
                end

                allowedEdges = [];

                % Горизонтальные связи
                for m = 1:M
                    base = (m-1)*N;
                    for k = 1:N
                        i = base + k;
                        j = base + mod(k,N) + 1;
                        allowedEdges = [allowedEdges; i j];
                    end
                end

                % Верхняя вершина соединяется с двумя нижними (лево/право)
                for m = 1:M-1
                    base     = (m-1)*N;
                    nextBase = m*N;
                    for k = 1:N
                        top   = nextBase + k;           % вершина верхнего уровня
                        left  = base + k;               % нижняя левая
                        right = base + mod(k,N) + 1;    % нижняя правая

                        allowedEdges = [allowedEdges;
                                        top left;
                                        top right];
                    end
                end

                % Удаляем дубликаты
                allowedEdges = unique(sort(allowedEdges,2), 'rows');

            else
                error('GridType должен быть "rect" или "tri".');
            end

            % --- Метаданные ---
            meta = struct( ...
                'ShapeType', 'cylinder', ...
                'PlacementType', 'surface', ...
                'GridType', gridType, ...
                'EdgeLength', L, ...
                'Description', templateParams.Description);
        end
    end

    methods (Static, Access = private)
        function N = estimateCylinderNSides(R, L)
            % Оценивает NSides по радиусу и желаемой длине ребра L.
            % Минимизируем |длина хорды - L|.

            if L <= 0 || R <= 0
                error('R и EdgeLength должны быть положительными.');
            end

            % Первичная оценка по длине дуги
            N_est = round(2*pi*R / L);
            N_est = max(N_est, 3);

            % Поиск в окрестности N_est
            candidates = max(3, N_est-5) : (N_est+5);
            candidates(candidates < 3) = [];

            errs = zeros(size(candidates));
            for i = 1:numel(candidates)
                Ntry = candidates(i);
                chord = 2*R*sin(pi/Ntry);
                errs(i) = abs(chord - L);
            end

            [~, idxMin] = min(errs);
            N = candidates(idxMin);
        end

        function M = estimateCylinderNLevels(H, L)
            % Оценивает NLevels по высоте H и желаемой длине ребра L.
            % Минимизируем |шаг по Z - L|.

            if L <= 0 || H <= 0
                error('Height и EdgeLength должны быть положительными.');
            end

            % Если один уровень — вообще без вертикальных рёбер (неинтересно),
            % поэтому минимум два уровня (один пояс вертикальных рёбер).
            M_est = round(H / L) + 1;
            M_est = max(M_est, 2);

            candidates = max(2, M_est-5) : (M_est+5);
            candidates(candidates < 2) = [];

            errs = zeros(size(candidates));
            for i = 1:numel(candidates)
                Mtry = candidates(i);
                dz = H / (Mtry - 1);
                errs(i) = abs(dz - L);
            end

            [~, idxMin] = min(errs);
            M = candidates(idxMin);
        end
    end
end
