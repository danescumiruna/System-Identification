function [MSE_Predictie, MSE_Simulare] = MSE_Validare(m_max, n_max, uid, yid, uval, yval)
    MSE_Predictie = zeros(m_max, n_max);
    MSE_Simulare = zeros(m_max, n_max);
    for m = 1:m_max
        for n = 1:n_max
            % Identificare
            clear phi;
            X = zeros(1, 2*n);
            for k = 1:length(uid)  
                % crearea vectorului de semnale intarziate
                col = 1;
                for i = 1:n
                    if k - i > 0
                        X(col) = yid(k - i);
                    else
                        X(col) = 0;
                    end
                    col = col + 1;
                end
                for i = 1:n
                    if k - i > 0
                        X(col) = uid(k - i);
                    else
                        X(col) = 0;
                    end
                    col = col + 1;
                end

                pows = zeros(1, length(X)); % initializarea vectorului de puteri
                phi(k, :) = back(1, pows, m, X); % generarea unei linii din matricea de regresori, apeland backtracking
                clear back; % folosit pentru stergerea variabilelor de tip persistent din functia back()
            end
            theta = phi\yid; % determinarea vectorului de parametri prin aplicarea regresiei liniare
            
            % Validare: Predictie
            X = zeros(1, 2*n);
            for k = 1:length(uval)
                % crearea vectorului de semnale intarziate
                col = 1;
                for i = 1:n
                    if k - i > 0
                        X(col) = yval(k - i);
                    else
                        X(col) = 0;
                    end
                    col = col + 1;
                end
                for i = 1:n
                    if k - i > 0
                        X(col) = uval(k - i);
                    else
                        X(col) = 0;
                    end
                    col = col + 1;
                end

                pows = zeros(1, length(X)); % initializarea vectorului de puteri
                phi(k, :) = back(1, pows, m, X); % generarea unei linii din matricea de regresori, apeland backtracking
                clear back; % folosit pentru stergerea variabilelor de tip persistent din functia back()
            end
            yhval = phi*theta; % determinarea semnalului de iesire aproximat folosind matricea de regresori si vectorul de parametri calculat pe datele de identificare
            MSE_Predictie(m, n) = sum((yhval - yval).^2)/length(uval); % calculare MSE
            
            % Validare: Simulare
            X = zeros(1, 2*n);
            yhval = zeros(length(uval), 1);
            for k = 1:length(uval)
                % crearea vectorului de semnale intarziate
                col = 1;
                for i = 1:n
                    if k - i > 0
                        X(col) = yhval(k - i);
                    else
                        X(col) = 0;
                    end
                    col = col + 1;
                end
                for i = 1:n
                    if k - i > 0
                        X(col) = uval(k - i);
                    else
                        X(col) = 0;
                    end
                    col = col + 1;
                end

                pows = zeros(1, length(X)); % initializarea vectorului de puteri
                phi(k, :) = back(1, pows, m, X); % generarea unei linii din matricea de regresori, apeland backtracking
                clear back; % folosit pentru stergerea variabilelor de tip persistent din functia back()

                yhval(k) = phi(k, :)*theta; % determinarea semnalului de iesire aproximat utilizand linia k a matricii de regresori si vectorul de parametri calculat in predictie pe datele de identificare
            end
            MSE_Simulare(m, n) = sum((yhval - yval).^2)/length(uval); % calculare MSE
        end
    end
end