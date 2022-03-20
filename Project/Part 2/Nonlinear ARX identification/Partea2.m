%% Date initiale
clear; clc; close all;
load('iddata-05.mat');
plot(id);
figure; plot(val);
%% Identificare: Predictie
clc; close all;
uid = id.U;
yid = id.Y;
na = 1;
nb = 1;
m = 5;

tid = 0:id.Ts:(length(uid) - 1)*id.Ts;
X = zeros(1, na + nb); 
for k = 1:length(uid) 
    % crearea vectorului de semnale intarziate
    col = 1;
    for i = 1:na
        if k - i > 0
            X(col) = yid(k - i);
        else
            X(col) = 0;
        end
        col = col + 1;
    end
    for i = 1:nb
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

MSE_id = sum((yhid - yid).^2)/length(uid); % calcularea MSE 
plot(tid, yid, tid, yhid,'Linewidth',2);
title(['MSE_{id} = ', num2str(MSE_id)], 'FontSize', 16);
xlabel('t_{id}'), ylabel('y');
legend('y_{id}','y_{id}^\^');
%% Identificare: Simulare
clc; close all;
tid = 0:id.Ts:(length(uid) - 1)*id.Ts;
X = zeros(1, na + nb);
yhid = zeros(length(uid), 1);
for k = 1:length(uid)
    % crearea vectorului de semnale intarziate
    col = 1;
    for i = 1:na
        if k - i > 0
            X(col) = yhid(k - i);
        else
            X(col) = 0;
        end
        col = col + 1;
    end
    for i = 1:nb
        if k - i > 0
            X(col) = uid(k - i);
        else
            X(col) = 0;
        end
        col = col + 1;
    end

    pows = zeros(1, length(X)); % initializarea vectorului de puteri
    phi(k, :) = back(1, pows, m, X);  % generarea unei linii din matricea de regresori, apeland backtracking
    clear back; % folosit pentru stergerea variabilelor de tip persistent din functia back()
    
    yhid(k) = phi(k, :)*theta; % determinarea semnalului de iesire aproximat utilizand linia k a matricii de regresori si vectorul de parametri calculat in predictie
end

MSE_id = sum((yhid - yid).^2)/length(uid); % calcularea MSE
plot(tid, yid, tid, yhid, 'Linewidth',2);
title(['MSE_{id} = ', num2str(MSE_id)],'FontSize', 16);
xlabel('t_{id}'), ylabel('y');
legend('y_{id}','y_{id}^\^');
%% Validare: Predictie
clear phi;
uval = val.U;
yval = val.Y;

tval = 0:val.Ts:(length(uval) - 1)*val.Ts;
X = zeros(1, na + nb);
for k = 1:length(uval)
    % crearea vectorului de semnale intarziate
    col = 1;
    for i = 1:na
        if k - i > 0
            X(col) = yval(k - i);
        else
            X(col) = 0;
        end
        col = col + 1;
    end
    for i = 1:nb
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

MSE_val = sum((yhval - yval).^2)/length(uval); % calcularea MSE
plot(tval, yval, tval, yhval, 'Linewidth',2);
title(['MSE_{val} = ', num2str(MSE_val)], 'FontSize', 16);
xlabel('t_{val}'), ylabel('y');
legend('y_{val}','y_{val}^\^');
%% Validare: Simulare
clc; close all;
X = zeros(1, na + nb);
yhval = zeros(length(uval), 1);
for k = 1:length(uval)
    % crearea vectorului de semnale intarziate
    col = 1;
    for i = 1:na
        if k - i > 0
            X(col) = yhval(k - i);
        else
            X(col) = 0;
        end
        col = col + 1;
    end
    for i = 1:nb
        if k - i > 0
            X(col) = uval(k - i);
        else
            X(col) = 0;
        end
        col = col + 1;
    end  
    
    pows = zeros(1, length(X)); % initializarea vectorului de puteri
    phi(k, :) = back(1, pows, m, X);% generarea unei linii din matricea de regresori, apeland backtracking
    clear back; % folosit pentru stergerea variabilelor de tip persistent din functia back()
    
    yhval(k) = phi(k, :)*theta; % determinarea semnalului de iesire aproximat utilizand linia k a matricii de regresori si vectorul de parametri calculat in predictie pe datele de identificare
end

MSE_val = sum((yhval - yval).^2)/length(uval); % calularea MSE
plot(tval, yval, tval, yhval, 'Linewidth',2);
title(['MSE_{val} = ', num2str(MSE_val)], 'FontSize', 16);
xlabel('t_{val}'), ylabel('y');
legend('y_{val}','y_{val}^\^');
%% MSE Validare
clc; close all;
m_max = 5; n_max = 3;
[MSE_Predictie, MSE_Simulare] = MSE_Validare(m_max, n_max, uid, yid, uval, yval);
%% MSE_min Validare
% determinarea ordinelor na, nb si a gradului m optime
m_min = 1; n_min = 1;
MSE_min = MSE_Predictie(m_min, n_min);
for m = 1:m_max
    for n = 1:n_max
        if MSE_Predictie(m, n) < MSE_min
            m_min = m;
            n_min = n;
            MSE_min = MSE_Predictie(m, n);
        end
    end
end

m_min = 1; n_min = 1;
MSE_min = MSE_Simulare(m_min, n_min);
for m = 1:m_max
    for n = 1:n_max
        if MSE_Simulare(m, n) < MSE_min
            m_min = m;
            n_min = n;
            MSE_min = MSE_Simulare(m, n);
        end
    end
end

%% MSE Identificare
clc; close all;
m_max = 5; n_max = 3;
[MSE_Predictie, MSE_Simulare] = MSE_Identificare(m_max, n_max, uid, yid);
%% MSE_min Identificare
% determinarea ordinelor na, nb si a gradului m optime
m_min = 1; n_min = 1;
MSE_min = MSE_Predictie(m_min, n_min);
for m = 1:m_max
    for n = 1:n_max
        if MSE_Predictie(m, n) < MSE_min
            m_min = m;
            n_min = n;
            MSE_min = MSE_Predictie(m, n);
        end
    end
end

m_min = 1; n_min = 1;
MSE_min = MSE_Simulare(m_min, n_min);
for m = 1:m_max
    for n = 1:n_max
        if MSE_Simulare(m, n) < MSE_min
            m_min = m;
            n_min = n;
            MSE_min = MSE_Simulare(m, n);
        end
    end
end
