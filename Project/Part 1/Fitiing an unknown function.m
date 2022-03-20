%% Identificare pentru m = 13
clear; clc; close all;
load('proj_fit_05.mat');
x1 = id.X{1};
x2 = id.X{2};
y = id.Y;
mesh(x1, x2, y);

m = 13;
nr = 0;

% Determinare numar de parametri
for a = 0:m
    for b = 0:m
        if(a+b <= m)
            nr = nr + 1;
        end
    end
end

% Construire matrice de regresori (phi) pentru datele de identificare
row = 1;
col = 1;
phi = zeros(id.dims(1)*id.dims(2), nr); % initializarea matricei phi
for i = 1:id.dims(1)
    for j = 1:id.dims(2) % dimensiunile matricei phi 
        for a = 0:m % a reprezinta exponentul coordonatei x1
            for b = 0:m  % b reprezinta exponentul coordonatei x2
                if(a+b <= m) % suma exponetilor este maxim m 
                    % fiecare element al matricei phi il calculam ca fiind
                    % produsul dintre cele doua coordonate ale intrarii la
                    % diferite puteri, acestia fiind regresorii
                    phi(row, col) =  x1(i)^a*x2(j)^b;
                    col = col + 1;
                end
            end
        end 
        col = 1;
        row = row+1;
    end
end
% Transformam matricea Y intr-o matrice coloana pentru a putea aplica regresia
yid = y(:);
theta = phi\yid; % aplicam regresia propriu-zisa
y_hatcol = phi*theta; % aproximam iesirea cu ajutorul matricii de regresori
                       % si a vectorului de parametri calculati anterior

% Redimensionam y_hatcol pentru a readuce matricea la forma initiala
y_hatId = zeros(id.dims(1), id.dims(2)); % initializarea matricii de aproximare
for i = 1:id.dims(1)
    % reconstruirea matricii din matrice coloana in matrice de
    % id.dims(1) x id.dims(2)
    y_hatId(:,i) = y_hatcol((i-1)*id.dims(1)+1:i*id.dims(1));
end

% Calculam MSE pe datele de identificare pentru m = 13
MSE = sum((y_hatcol - yid).^2)/length(yid);

% Suprapunem modelul gasit pe datele de identificare cu sistemul real
hold on
mesh(x1, x2, y_hatId);
title(['MSE = ', num2str(MSE)],'FontSize',16);
xlabel('x_1'); ylabel('x_2'); zlabel('y');
hold off
%% Validare pentru m = 13
clc; close all;
x1 = val.X{1};
x2 = val.X{2};
y = val.Y;
mesh(x1, x2, y);

% Construire matrice de regresori (phi) pentru datele de validare
row = 1;
col = 1;  
phi = zeros(val.dims(1)*val.dims(2), nr); % initializarea matricei phi

for i = 1:val.dims(1)
    for j = 1:val.dims(2) % dimensiunile matricei phi
        for a = 0:m % a reprezinta exponentul coordonatei x1
            for b = 0:m % b reprezinta exponentul coordonatei x2
                if(a+b <= m) % suma exponetilor este maxim m
                    % fiecare element al matricei phi il calculam ca fiind
                    % produsul dintre cele doua coordonate ale intrarii la
                    % diferite puteri, acestia fiind regresorii
                    phi(row, col) =  x1(i)^a*x2(j)^b;
                    col = col + 1;
                end
            end
        end 
        col = 1;
        row = row+1;
    end
end

% Transformam matricea Y intr-o matrice coloana pentru a putea aplica regresia
yval = y(:);
y_hatVal = zeros(val.dims(1), val.dims(2)); %initializarea matricii de aproximare

% Folosim matricea de parametri (theta) determinata anterior pe datele de
   % identificare
y_hatcol = phi*theta;  % aproximam iesirea cu ajutorul matricii de regresori
                       % calculata pe datele de validare si a vectorului de
                       % parametri calculat pe datele de identificare

% Redimensionam y_hatcol pentru a readuce matricea la forma initiala
for i = 1:val.dims(1)
    % reconstruirea matricii din matrice coloana in matrice de
    % val.dims(1) x val.dims(2)
    y_hatVal(:,i) = y_hatcol((i-1)*val.dims(1)+1:i*val.dims(1));
end
% Calculam MSE pe datele de validare pentru m = 13
MSE = sum((y_hatcol - yval).^2)/length(yval);

% Suprapunem modelul gasit pe datele de validare cu sistemul real
hold on
mesh(x1, x2, y_hatVal);
title(['MSE = ', num2str(MSE)],'FontSize',16)
xlabel('x_1'); ylabel('x_2'); zlabel('y');
hold off
%% MSE in functie de m 
clc; clear; close all;
load('proj_fit_05.mat');

for m = 1:25
 % Determinare numar de parametri
    nr = 0;
    for a = 0:m
        for b = 0:m
            if(a+b <= m)
                nr = nr + 1;
            end
        end
    end
    
 % Construire matrice de regresori (phi) pentru datele de identificare
    row = 1;
    col = 1;
    phi = zeros(id.dims(1)*id.dims(2), nr); % initializarea matricei phi
    x1 = id.X{1};
    x2 = id.X{2};
    y = id.Y;
    for i = 1:id.dims(1)
        for j = 1:id.dims(2) % dimensiunile matricei phi 
            for a = 0:m % a reprezinta exponentul coordonatei x1
                for b = 0:m % b reprezinta exponentul coordonatei x2
                    if(a+b <= m) % suma exponetilor este maxim m
                        % fiecare element al matricei phi il calculam ca fiind
                        % produsul dintre cele doua coordonate ale intrarii la
                        % diferite puteri, acestia fiind regresorii
                        phi(row, col) =  x1(i)^a*x2(j)^b;
                        col = col + 1;
                    end
                end
            end 
            col = 1;
            row = row+1;
        end
    end
 % Transformam matricea Y intr-o matrice coloana pentru a putea aplica regresia
    yid = y(:);
    theta = phi\yid;  % aplicam regresia propriu-zisa
    y_hatcol = phi*theta; % aproximam iesirea cu ajutorul matricii de regresori
                            % si a vectorului de parametri calculati anterior
    
    y_hatId = zeros(id.dims(1), id.dims(2)); % initializarea matricii de aproximare
    for i = 1:id.dims(1)
        % reconstruirea matricii din matrice coloana in matrice de
        % id.dims(1) x id.dims(2)
        y_hatId(:,i) = y_hatcol((i-1)*id.dims(1)+1:i*id.dims(1));
    end
    
    % Calculam MSE pe datele de identificare pentru m variabil
    MSEid(m) = sum((y_hatcol - yid).^2)/length(yid);
    
 % Construire matrice de regresori (phi) pentru datele de validare   
    row = 1;
    col = 1;
    phi = zeros(val.dims(1)*val.dims(2), nr); % initializarea matricei phi
    x1 = val.X{1};
    x2 = val.X{2};
    y = val.Y;
    
    for i = 1:val.dims(1)
        for j = 1:val.dims(2) % dimensiunile matricei phi
            for a = 0:m % a reprezinta exponentul coordonatei x1
                for b = 0:m % b reprezinta exponentul coordonatei x2
                    if(a+b <= m) % suma exponetilor este maxim m
                        % fiecare element al matricei phi il calculam ca fiind
                        % produsul dintre cele doua coordonate ale intrarii la
                        % diferite puteri, acestia fiind regresorii
                        phi(row, col) =  x1(i)^a*x2(j)^b;
                        col = col + 1;
                    end
                end
            end 
            col = 1;
            row = row+1;
        end
    end
 % Transformam matricea Y intr-o matrice coloana pentru a putea aplica regresia
    yval = y(:);
    y_hatcol = phi*theta; % aproximam iesirea cu ajutorul matricii de regresori
                          % calculata pe datele de validare si a vectorului de
                          % parametri calculat pe datele de identificare
    y_hatVal = zeros(val.dims(1), val.dims(2)); %initializarea matricii de aproximare
    
% Redimensionam matricea y_hatcol pentru a o readuce la forma initiala
    for i = 1:val.dims(1)
        % reconstruirea matricii din matrice coloana in matrice de
        % val.dims(1) x val.dims(2)
        y_hatVal(:,i) = y_hatcol((i-1)*val.dims(1)+1:i*val.dims(1));
    end
    
 % Calculam MSE pe datele de validare pentru m din intervalul [1,25]
    MSEval(m) = sum((y_hatcol - yval).^2)/length(yval);
end

% Afisam graficul lui MSE din intervalul [1,25] atat pentru datele de
% identificare, cat si pentru datele de validare
figure, plot(MSEid,'LineWidth',2);
xlabel('m'); ylabel('MSE');
title('MSE identificare','FontSize',16); 
figure, plot(MSEval,'LineWidth',2);
title('MSE validare','FontSize',16);
xlabel('m'); ylabel('MSE');
MSEmin = MSEval(1);
mmin = 1;
for i = 2:length(MSEval) % calculul gradului minim
    if(MSEval(i) < MSEmin)
        MSEmin = MSEval(i);
        mmin = i;
    end
end

hold on
plot(mmin, MSEmin, 'r*','LineWidth',2); % evidentierea gradului minim 
hold off