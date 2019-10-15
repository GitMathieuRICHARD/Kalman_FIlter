%-----------------------------------------------------------
%-------------------------- KALMAN Filter ------------------
%-----------------------------------------------------------
%-----------------------------------------------------------
%------ Cette algorithme a pour fonction d'exposer le ------
%------ fonctionnement du filtre de Kalman à état ----------
%------ stationnaire. Il peu estimer la sortie y en --------
%------ fonction des mesures bruitées. ---------------------
%-----------------------------------------------------------
%------- Cette fonction détermine le gain optimal du -------
%------- filtre en régime permanent M en fonction de la ----
%------- covariance du bruit de processus Q et de la -------
%------- covariance de bruit du capteur R.------------------
%-----------------------------------------------------------
%-----------------------------------------------------------
%----------------- Edité par Mathieu RICHARD ---------------
%-----------------------------------------------------------
%------------------------------- Le 15/10/2019 -------------
%-----------------------------------------------------------

% Lien 1 : https://fr.mathworks.com/help/control/examples/kalman-filter-design.html
% Lien 2 : https://fr.mathworks.com/help/control/ug/kalman-filtering.html

%% Steady-State Kalman Filter Design

%Test Kalman Filter

A = [1.1269   -0.4940    0.1129,
     1.0000         0         0,
          0    1.0000         0];

B = [-0.3832
      0.5919
      0.5191];

C = [1 0 0];

D = 0;

% Spécifiez d'abord le modèle plant + bruit. ATTENTION: réglez la durée d'échantillonnage sur -1 pour marquer l'installation comme étant discrète.
% La fonction ss permet de
Plant = ss(A,[B B],C,0,-1,'inputname',{'u' 'w'},'outputname','y')

Q = 2.3; % Un nombre plus grand que 0

R = 1; % Un nombre plus grand que 0

[kalmf,L,~,M,Z] = kalman(Plant,Q,R);

kalmf = kalmf(1,:);

M,   % innovation gain

% Premièrement, il faut construire un modele de plant complète avec u,w,v en entrée et y et yv en sortie :
a = A;
b = [B B 0*B];
c = [C;C];
d = [0 0 0;0 0 1];
P = ss(a,b,c,d,-1,'inputname',{'u' 'w' 'v'},'outputname',{'y' 'yv'});

sys = parallel(P,kalmf,1,1,[],[]);

SimModel = feedback(sys,1,4,2,1);
SimModel = SimModel([1 3],[1 2 3]);     % Delete yv form I/O

SimModel.inputname

SimModel.outputname

t = (0:100)';
u = sin(t/5);

rng(10,'twister');
w = sqrt(Q)*randn(length(t),1);
v = sqrt(R)*randn(length(t),1);

out = lsim(SimModel,[w,v,u]);

y = out(:,1);   % true response
ye = out(:,2);  % filtered response
yv = y + v;     % measured response

clf
subplot(211), plot(t,y,'b',t,ye,'r--'),
xlabel('No. of samples'), ylabel('Output')
title('Kalman filter response')
subplot(212), plot(t,y-yv,'g',t,y-ye,'r--'),
xlabel('No. of samples'), ylabel('Error')

MeasErr = y-yv;
MeasErrCov = sum(MeasErr.*MeasErr)/length(MeasErr);
EstErr = y-ye;
EstErrCov = sum(EstErr.*EstErr)/length(EstErr);

MeasErrCov

EstErrCov

%% Time-Varying Kalman Filter Design

sys = ss(A,B,C,D,-1);
y = lsim(sys,u+w);   % w = process noise
yv = y + v;          % v = meas. noise

P=B*Q*B';         % Initial error covariance
x=zeros(3,1);     % Initial condition on the state
ye = zeros(length(t),1);
ycov = zeros(length(t),1);
errcov = zeros(length(t),1);

for i=1:length(t)
  % Measurement update
  Mn = P*C'/(C*P*C'+R);
  x = x + Mn*(yv(i)-C*x);  % x[n|n]
  P = (eye(3)-Mn*C)*P;     % P[n|n]

  ye(i) = C*x;
  errcov(i) = C*P*C';

  % Time update
  x = A*x + B*u(i);        % x[n+1|n]
  P = A*P*A' + B*Q*B';     % P[n+1|n]
end

subplot(211), plot(t,y,'b',t,ye,'r--'),
xlabel('No. of samples'), ylabel('Output')
title('Response with time-varying Kalman filter')
subplot(212), plot(t,y-yv,'g',t,y-ye,'r--'),
xlabel('No. of samples'), ylabel('Error')

subplot(211)
plot(t,errcov), ylabel('Error Covar'),

MeasErr = y-yv;
MeasErrCov = sum(MeasErr.*MeasErr)/length(MeasErr);
EstErr = y-ye;
EstErrCov = sum(EstErr.*EstErr)/length(EstErr);

MeasErrCov

EstErrCov

M,Mn




