function [finalF] = modify_f(F)
% Projekt 1, zadanie 24
% Łukasz Szymczyk, 320744
%
% Funkcja zamienia wektor funkcji anonimowych zmiennej x postaci 
% {b, a_0, ... , a_m} na funkcję dwóch zmiennych x,y.
% Wejście:
%   F     - wektor uchwytów do funkcji zmiennej x w postaci
%           {b, a_0, ... , a_m}
% Wyjście:
%   finalF - funkcja anonimowa postaci finalF = @(x,y) 
%            [y(2);...;y(m);f(x,y)]
%            gdzie x to pojedyncza wartość, y - wektor, a
%            f(x,y) = (b(x) - a_0(x)*y(1) - ... - a_m-1(x)*y(m))/a_m(x)

% Zainicjowanie potrzebnych zmiennych
N = size(F,2);
f = F{N};
Z = @(x,y) F{1}(x);

% Wyznaczenie f(x,y)
for i = 2:N-1
    G = @(x,y) F{i}(x)*y(i-1);
    Z = @(x,y) Z(x,y)-G(x,y);
end
Z = @(x,y) Z(x,y)/f(x);

% Przypadek, gdy równanie jest 1-go rzędu
if (N == 3)
    finalF = Z;
    return
end

% Tworzenie funkcji [y(2);...;y(m)]
H = @(x,y) y(2);
for i = 2:N-3
    H = @(x,y) [H(x,y);y(i+1)];
end

% Stworzenie [y(2);...;y(m);f(x,y)]
finalF = @(x,y) [H(x,y);Z(x,y)];

end % function

