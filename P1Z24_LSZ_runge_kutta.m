function [X,Y,h] = P1Z24_LSZ_runge_kutta(F,range,Y0,N)
% Projekt 1, zadanie 24
% Łukasz Szymczyk, 320744
%
% Rozwiązywanie liniowych równań różniczkowych dowolnego rzędu,
% postaci a_m(x)*y^(m) + a_m-1(x)*y(m-1)+...+a_1(x)*y'+a_0(x)*y = b(x),
% iteracyjną metodą Rungego-Kutty rzędu 4-go (wzór 3/8)
% Wejście:
%   F     - wektor uchwytów do funkcji zmiennej x w postaci
%           {b, a_0, ... , a_m}
%   range - dwuelementowy wektor [a, b] określający przedział, na którym
%           przybliżamy szukaną funkcję
%   Y0    - m-elementowy wektor zawierający warunki początkowe równania
%           [y0, y'0, ... , ]
%   N     - liczba kroków
% Wyjście:
%   X     - wektor (rozmiar N+1x1) równoodległych węzłów, w których 
%           przybliżana jest wartość szukanej funkcji
%   Y     - wektor (rozmiar N+1x1) zawierający wyznacznone przybliżenia 
%           wartości szukanej funkcji w węzłach zawartych w wektorze X
%   h     - krok całkowania

f = modify_f(F);
[X,Y,h] = runge_kutta_3_8(f,range,Y0,N);

end % function

