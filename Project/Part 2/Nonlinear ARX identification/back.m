function phik = back(k, pows, m, X)
    % persistent face ca variabilele sa si pastreze valoarea de la un apel
    % la altul
    persistent c
    persistent lin
    if (isempty(c)) % initializare c
        c = 1;
    end
    if k > length(pows)
        p = 1;
        for i = 1:length(X)
            p = p*X(i)^pows(i); % calcularea produsului elementelor vectorului X la diferite puteri
        end 
        lin(c) = p;
        c = c + 1;
        phik = lin; % construirea liniei k a matricii de regresori
    else
        for i = 0:m
            pows(k) = i;
            s = 0;
            for j = 1:k
                s = s + pows(j); % calculare grad curent
            end
            if s <= m % verificare grad curent < m 
              phik = back(k + 1, pows, m, X); % se apeleaza functia recursiv pentru pasul urmator
            end
            pows(k) = 0;
        end
    end
end