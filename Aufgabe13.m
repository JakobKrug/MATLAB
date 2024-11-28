close all; % Schließt alle offenen Figurenfenster
clear vars; % Löscht alle Variablen aus dem Workspace

%% 13d)
n = 10; % Setzt die Größe der Matrix auf 10
diagonals = {rand(n - 1, 1), rand(n, 1), rand(n - 1, 1)}; % Erstellt Zufallswerte für die Unter-, Haupt- und Oberdiagonale
b = rand(n, 1); % Erstellt einen Zufallsvektor b

% Testet die LR-Zerlegung mit den generierten Diagonalen und dem Vektor b
test_lr_tridiag(diagonals{1}, diagonals{2}, diagonals{3}, b);

% Setzt den dritten Wert der Hauptdiagonale auf 0, um einen speziellen Fall zu testen
diagonals{2}(3) = 0;
% Testet die LR-Zerlegung erneut mit der modifizierten Hauptdiagonale und dem Vektor b
test_lr_tridiag(diagonals{1}, diagonals{2}, diagonals{3}, b);

% Funktion zur Durchführung der LR-Zerlegung einer Tridiagonalmatrix
function [lund, rhd, rond] = my_lr_tridiag(und, hd, ond)
    n = length(hd); % Bestimmt die Länge der Hauptdiagonale
    lund = zeros(n - 1, 1); % Initialisiert die Unterdiagonale der L-Matrix
    rhd = hd; % Kopiert die Hauptdiagonale in die R-Matrix
    rond = ond; % Kopiert die Oberdiagonale in die R-Matrix

    % Überprüft, ob es eine Null auf der Hauptdiagonale gibt
    if any(rhd == 0)
        error('Error: 0 on main diag encountered, cannot perform LU-decomposition w/o pivoting'); % Gibt einen Fehler aus, wenn eine Null gefunden wird
    end

    % Führt die LR-Zerlegung durch
    for i = 1:n - 1
        lund(i) = und(i) / rhd(i); % Berechnet den Wert der Unterdiagonale der L-Matrix
        rhd(i + 1) = rhd(i + 1) - lund(i) * rond(i); % Aktualisiert die Hauptdiagonale der R-Matrix
    end

end

% Funktion zum Testen der LR-Zerlegung einer Tridiagonalmatrix
function test_lr_tridiag(und, hd, ond, b)

    try
        % Führt die LR-Zerlegung durch und gibt die Matrizen aus
        [lund, rhd, rond] = my_lr_tridiag(und, hd, ond);
        disp('Decomposition of:');
        disp(diag(und, -1) + diag(hd) + diag(ond, 1)); % Zeigt die ursprüngliche Matrix
        disp(diag(lund, -1) + diag(rhd) + diag(rond, 1)); % Zeigt die zerlegte Matrix

        % Löst das Gleichungssystem und zeigt die Lösung
        x = solve_lr_tridiag(lund, rhd, rond, b);
        A = diag(lund, -1) + diag(rhd) + diag(rond, 1); % Rekonstruiert die Matrix aus den Diagonalen
        disp(x); % Zeigt die berechnete Lösung
        disp(A \ b); % Zeigt die Lösung durch direkte Berechnung
    catch e
        disp(e.message); % Gibt eine Fehlermeldung aus, falls ein Fehler auftritt
    end

end

% Funktion zum Lösen eines Gleichungssystems mit einer LR-Zerlegung
function x = solve_lr_tridiag(lund, rhd, rond, b)
    n = size(rhd, 1); % Bestimmt die Größe des Systems

    % L*y=b
    y = zeros(n, 1);
    y(1) = b(1); % Initialisiert die erste Komponente von y

    for i = 2:n
        y(i) = b(i) - lund(i - 1) * y(i - 1); % Führt die Vorwärtssubstitution durch
    end

    % R*x=y
    x = zeros(n, 1);
    x(n) = y(n) / rhd(n); % Initialisiert die letzte Komponente von x

    for i = n - 1:-1:1
        x(i) = (y(i) - x(i + 1) * rond(i)) / rhd(i); % Führt die Rückwärtssubstitution durch
    end

end
