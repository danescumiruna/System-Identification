function [MSE_Predictie, MSE_Simulare] = MSE_Identificare(m_max, n_max, uid, yid)
    MSE_Predictie = zeros(m_max, n_max);
    MSE_Simulare = zeros(m_max, n_max);
    for m = 1:m_max
        for n = 1:n_max
            % Identificare: Predictie
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
            yhid = phi*theta; % determinarea semnalului de iesire aproximat
            MSE_Predictie(m, n) = sum((yhid - yid).^2)/length(uid); % calculare MSE
            
            % Identificare: Simulare
            X = zeros(1, 2*n);
            yhid = zeros(length(uid), 1);
            for k = 1:length(uid)    
                % crearea vectorului de semnale intarziate
                col = 1;
                for i = 1:n
                    if k - i > 0
                        X(col) = yhid(k - i);
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

                yhid(k) = phi(k, :)*theta;  % determinarea semnalului de iesire aproximat utilizand linia k a matricii de regresori si vectorul de parametri calculat in predictie
            end
            MSE_Simulare(m, n) = sum((yhid - yid).^2)/length(uid); % calculare MSE
        end
    end
end