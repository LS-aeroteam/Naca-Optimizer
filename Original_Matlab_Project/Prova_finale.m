clear all; close all; clc;

%% === Parametri utente ===
naca      = '0012';     % Esempio: '2412'
alpha     =  0;         % Incidenza geometrica(gradi)
nPanels   = 160;        % Numero di pannelli totali 

%% ===Impostazioni codice===
plotStreamlines = true;      % true = calcolo streamlines, false = skip per esecuzione veloce
plotCp= true;                % true = calcolo streamlines, false = skip per esecuzione veloce



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Generazione Profilo

% Parametri NACA 4 cifre 
M = str2double(naca(1))/100;    % Max camber
P = str2double(naca(2))/10;     % Posizione max camber
t = str2double(naca(3:4))/100;  % Spessore relativo

% Coordinate con distribuzione cosinusoidale 
nHalf = nPanels/2;                              % punti per semiprofilo
beta  = linspace(0, pi, nHalf+1);               
x     = (1 - cos(beta))/2 ;                     % distribuzione cosinusoidale

% Spessore profilo
yt = 5*t*( 0.2969*sqrt(x) - 0.1260*(x) ...
           - 0.3516*(x).^2 + 0.2843*(x).^3 ...
           - 0.1015*(x).^4 );

% Linea media (camber) e sua derivata
yc  = zeros(size(x));
dyc = zeros(size(x));

if M ~= 0
    for i = 1:length(x)
        xi = x(i);
        if xi < P
            yc(i)  = (M/P^2) * (2*P*xi - xi^2);
            dyc(i) = (2*M/P^2) * (P - xi);
        else
            yc(i)  = M/(1-P)^2 * ((1-2*P) + 2*P*xi - xi^2);
            dyc(i) = 2*M/(1-P)^2 * (P - xi);
        end
    end
end

theta = atan(dyc);    % angolo linea media con orizzonte

% Punti del profilo dorso (upper) e ventre (lower)
xu = x - yt.*sin(theta);
yu = yc + yt.*cos(theta);
xl = x + yt.*sin(theta);
yl = yc - yt.*cos(theta);

%% Costruzione pannelli e checkpoints

% Punti in senso orario: ventre (BU→BA) + dorso (BA→BU)
XB = [xl(end:-1:1), xu(2:end)];
YB = [yl(end:-1:1), yu(2:end)];

% Calcolo punti di controllo, lunghezze e angoli con orizzonte dei pannelli
nPanels = length(XB)-1;           % Aggiornato per sicurezza
XC = zeros(1, nPanels);
YC = zeros(1, nPanels);
S  = zeros(1, nPanels);
PSI = zeros(1, nPanels);

for i = 1:nPanels
    XC(i)  = (XB(i) + XB(i+1))/2;
    YC(i)  = (YB(i) + YB(i+1))/2;
    S(i)   = hypot(XB(i+1)-XB(i), YB(i+1)-YB(i));
    PSI(i) = atan2d(YB(i+1)-YB(i), XB(i+1)-XB(i));
    if PSI(i) < 0, PSI(i) = PSI(i) + 360; end
end





%%  Matrici influenza per Vortex

K = zeros(nPanels);                                                   % Matrice velocità normali
L = zeros(nPanels);                                                   % Matrice velocità tangenziali

for i = 1:nPanels
    for j = 1:nPanels
        if j ~= i

            % Calcolo termini risoluzione analitica
            A  = -(XC(i)-XB(j))*cosd(PSI(j)) - (YC(i)-YB(j))*sind(PSI(j));
            B  = (XC(i)-XB(j))^2 + (YC(i)-YB(j))^2;
            Cn = -cosd(PSI(i)-PSI(j));
            Dn = (XC(i)-XB(j))*cosd(PSI(i)) + (YC(i)-YB(j))*sind(PSI(i));
            Ct = sind(PSI(j)-PSI(i));                                         
            Dt = (XC(i)-XB(j))*sind(PSI(i))-(YC(i)-YB(j))*cosd(PSI(i));       
            E  = sqrt(max(B-A^2,0));

            % Calcolo K(i,j)
            term1 = 0.5*Cn*log((S(j)^2+2*A*S(j)+B)/B);
            term2 = ((Dn-A*Cn)/E) * (atan2(S(j)+A, E) - atan2(A, E));
            K(i,j) = term1 + term2;

            % Calcolo L(i,j)
            term1  = 0.5*Ct*log((S(j)^2+2*A*S(j)+B)/B);                     
            term2  = ((Dt-A*Ct)/E)*(atan2((S(j)+A),E)-atan2(A,E));          
            L(i,j) = term1 + term2;                                         
        end
    end
end

%%  Matrici influenza per Source
I = zeros(nPanels);                                                   % Matrice velocità normali
J = zeros(nPanels);                                                   % Matrice velocità tangenziali


for i = 1:nPanels                                                          
    for j = 1:nPanels                                                      
        if (j ~= i)                                                         
            % Termini risoluzione analitica
            A  = -(XC(i)-XB(j))*cosd(PSI(j))-(YC(i)-YB(j))*sind(PSI(j));      
            B  = (XC(i)-XB(j))^2+(YC(i)-YB(j))^2;                           
            Cn = sind(PSI(i)-PSI(j));                                        
            Dn = -(XC(i)-XB(j))*sind(PSI(i))+(YC(i)-YB(j))*cosd(PSI(i));      
            Ct = -cosd(PSI(i)-PSI(j));                                       
            Dt = (XC(i)-XB(j))*cosd(PSI(i))+(YC(i)-YB(j))*sind(PSI(i));       
            E  = sqrt(max(B-A^2,0));                                              
            
            
            % Calcolo I(i,j) 
            term1  = 0.5*Cn*log((S(j)^2+2*A*S(j)+B)/B);                     
            term2  = ((Dn-A*Cn)/E)*(atan2((S(j)+A),E) - atan2(A,E));        
            I(i,j) = term1 + term2;                                         
            
            % Calcolo J(i,j) 
            term1  = 0.5*Ct*log((S(j)^2+2*A*S(j)+B)/B);                     
            term2  = ((Dt-A*Ct)/E)*(atan2((S(j)+A),E) - atan2(A,E));        
            J(i,j) = term1 + term2;                                         
        end
    end
end



%% Matrice influenza finale

A = zeros(nPanels,nPanels);                                                   
for i = 1:nPanels                                                          
    for j = 1:nPanels                                                     
        if (j == i)                                                         
            A(i,j) = pi;     % gestione "auto-influenza"                                                  
        else                                                                
            A(i,j) = I(i,j);                                                
        end
    end
end

% Colonna a destra per il gamma
for i = 1:nPanels                                                         
    A(i,nPanels+1) = sum(K(i,:));                                           
end

% Implementazione condizione di Kutta
for j = 1:1:nPanels                                                          
    A(nPanels+1,j) = (J(1,j) + J(nPanels,j));                                 % Contributo Source condizione Kutta 
end
A(nPanels+1,nPanels+1) = sum(L(1,:) + L(nPanels,:)) - 2*pi;                   % Contributo Vortex condizione Kutta 


% Vettore termini noti (flusso asintotico)
b = zeros(nPanels,1);                                                        
for i = 1:nPanels                                                          
    b(i) = -2*pi*sind(alpha-PSI(i));                                         
end

% Aggiungo ultimo elemento di b (condizione Kutta) 
b(nPanels+1) = -2*pi*(cosd(alpha-PSI(1)) + cosd(alpha-PSI(nPanels)));                

% Ridoluzione sistema
resArr = A\b;                                                               

% Separazione dei lambda dal gamma 
lambda = resArr(1:end-1);                                                   
gamma  = resArr(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calcolo Cl

Perimeter=sum(S);
GAMMA=gamma*Perimeter;
CL=-2*GAMMA;

%%  Risultati intermedi

fprintf('\n=== Recap ===\n');
fprintf('Profilo NACA: %s\n', naca);
fprintf('Incidenza geometrica: %.2f gradi\n', alpha);
fprintf('Numero pannelli: %d\n', nPanels)
fprintf('CL calcolato: %.4f\n', CL);


%% Calcolo velocità tangenziali e Cp sui pannelli
Vt = zeros(nPanels,1);                                                       
Cp = zeros(nPanels,1);                                                       
for i = 1:1:nPanels
    term1 = cosd(alpha-PSI(i));                                              % Termine flusso uniforme
    term2 = (1/(2*pi))*sum(lambda.*J(i,:)');                                % Termine delle sorgenti
    term3 = -gamma/2;                                                        % termine auto-influenza vortici
    term4 = (gamma/(2*pi))*sum(L(i,:));                                    % termine vortici
    
    Vt(i) = term1 + term2 + term3 + term4;                                  
    Cp(i) = 1-((Vt(i))^2);                                               
end






%% Disegno streamlines
if plotStreamlines==true 

    % Parametri griglia
    nGridX = 100;                                                           
    nGridY = 100;                                                           
    xVals  = [min(XB)-0.5 max(XB)+0.5];                                     %  [min, max] griglia x
    yVals  = [min(YB)-0.3 max(YB)+0.3];                                     %  [min, max] griglia y
    
    % Parametri streamlines
    stepsize = 0.01;                                                       
    slPct    = 25;                                                          % Percentage of streamlines of the grid
    Ysl      = linspace(yVals(1),yVals(2),floor((slPct/100)*nGridY))';      % Create array of Y streamline starting points
    
    % Generate the grid points
    Xgrid   = linspace(xVals(1),xVals(2),nGridX)';                          % griglia relativa ai valori x
    Ygrid   = linspace(yVals(1),yVals(2),nGridY)';                          % griglia realtiva ai valori y
    [XX,YY] = meshgrid(Xgrid,Ygrid);                                        
    
    % Campo di velocità
    Vx = zeros(nGridX,nGridY);                                              
    Vy = zeros(nGridX,nGridY);                                              
    


    % Generazione punti P=(XP;YP)
    for m = 1:1:nGridX
        for n = 1:1:nGridY
            XP      = XX(m,n);                                              % x del generico punto P
            YP      = YY(m,n);                                              % y del generico punto P
            

% Calcolo matrici influenza per streamlines (Source)
Mx = zeros(nPanels,1);                                                       
My = zeros(nPanels,1);                                                       


for j = 1:1:nPanels                                                         
    % Termini risoluzione analitica
    A  = -(XP-XB(j))*cosd(PSI(j))-(YP-YB(j))*sind(PSI(j));                   
    B  = (XP-XB(j))^2+(YP-YB(j))^2;                                         
    Cx = -cosd(PSI(j));                                                      
    Dx = XP - XB(j);                                                        
    Cy = -sind(PSI(j));                                                      
    Dy = YP - YB(j);                                                        
    E  = sqrt(max(0,(B-A^2)));                                                       
    
    
    
    term1 = 0.5*Cx*log((S(j)^2+2*A*S(j)+B)/B);                              
    term2 = ((Dx-A*Cx)/E)*(atan2((S(j)+A),E) - atan2(A,E));                 
    Mx(j) = term1 + term2;                                                  
    
    
    term1 = 0.5*Cy*log((S(j)^2+2*A*S(j)+B)/B);                              
    term2 = ((Dy-A*Cy)/E)*(atan2((S(j)+A),E) - atan2(A,E));                 
    My(j) = term1 + term2;                                                  
    
    
    
end

% Calcolo matrici influenza per streamlines (Vortex)
Nx = zeros(nPanels,1);                                                       
Ny = zeros(nPanels,1);                                                       


for j = 1:1:nPanels                                                          
   % Termini risoluzione analitica
    A  = -(XP-XB(j))*cosd(PSI(j))-(YP-YB(j))*sind(PSI(j));                    
    B  = (XP-XB(j))^2+(YP-YB(j))^2;                                         
    Cx = sind(PSI(j));                                                       
    Dx = -(YP-YB(j));                                                       
    Cy = -cosd(PSI(j));                                                      
    Dy = XP-XB(j);                                                          
    E  = sqrt(max(0,(B-A^2)));                                                       
    
   
    term1 = 0.5*Cx*log((S(j)^2+2*A*S(j)+B)/B);                              
    term2 = ((Dx-A*Cx)/E)*(atan2((S(j)+A),E) - atan2(A,E));                 
    Nx(j) = term1 + term2;                                                  
   
    term1 = 0.5*Cy*log((S(j)^2+2*A*S(j)+B)/B);                              
    term2 = ((Dy-A*Cy)/E)*(atan2((S(j)+A),E) - atan2(A,E));                 
    Ny(j) = term1 + term2;                                                  
    
	
end


% Velocità nulle dentro al profilo            
            [in,on] = inpolygon(XP,YP,XB,YB);
             if (in == 1 || on == 1)                                         % Se il punto è fuori dal poligono
                 Vx(m,n) = 0;                                                
                 Vy(m,n) = 0;                                                
             else                                                            
                Vx(m,n) = cosd(alpha) + sum(lambda.*Mx./(2*pi)) + ...    
                            sum(gamma.*Nx./(2*pi));
                Vy(m,n) = sind(alpha) + sum(lambda.*My./(2*pi)) + ...    
                            sum(gamma.*Ny./(2*pi));
            end
        end
     end
    
   
   


%% Plot Streamlines 
figure;
hold on;
axis equal;
xlabel('x');
ylabel('y');
title('Streamlines attorno al profilo','FontSize',18,'FontWeight','bold');

% Disegno del profilo
 patch(XB, YB, [1 1 1], 'EdgeColor', 'k', 'LineWidth', 1.5);

% Punti di partenza per le streamlines (verticale a sinistra)
startX = xVals(1) * ones(size(Ysl));
startY = Ysl;

% Disegno delle streamlines
streamline(XX, YY, Vx, Vy, startX, startY);

% Limiti grafico
xlim(xVals');
ylim(yVals');
grid on;

end  



%%  Plot Cp su dorso e ventre 
if plotCp==true


% Suddivisione pannelli appartenenti al ventre e dorso 
nHalf = nPanels/2;
xLower = XC(1:nHalf);
CpLower = Cp(1:nHalf);

xUpper = XC(nHalf+1:end);
CpUpper = Cp(nHalf+1:end);

% Grafico Cp
figure;
hold on; grid on;
plot(xUpper, CpUpper, '-', 'LineWidth', 1.5, 'DisplayName', 'Dorso');
plot(xLower, CpLower, '-', 'LineWidth', 1.5, 'DisplayName', 'Ventre');


set(gca, 'YDir','reverse'); 
xlabel('x/c');
ylabel('C_p');
title(['Distribuzione C_p attorno al NACA ', naca, ', \alpha = ', num2str(alpha), '°']);
legend('Location','best');

end