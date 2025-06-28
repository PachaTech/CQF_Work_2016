% Euler and Milstein discretization for Black-Scholes.

% Inline function for the Black Scholes call
C = inline('s*exp(-q*T)*normcdf((log(s/K) + (r-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-r*T)*normcdf((log(s/K) + (r-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T))',...
	's','K','r','q','v','T');

% Inline function for the Black Scholes put
P = inline('K*exp(-r*T)*normcdf(-(log(s/K) + (r-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (r-q+v^2/2)*T)/v/sqrt(T))',...
	's','K','r','q','v','T');
%===========================================
% Define the Option Input variables
%===========================================
S0        = 100;               % Spot price
K         = 100;               % Strike price
r         = 0.05;              % Risk free rate
q         = 0.0;               % Dividend yield
v         = 0.20;              % Volatility
mat       = 1.0;               % Time to maturity in years
PutCall   = 'C';               % 'P'ut or 'C'all
Averaging = 'A';               % 'A'rithmetic or 'G'eometric

%===========================================
% Monte Carlo Simulation settings.
%===========================================

N = 1000;             % Number of simulations.
T = 1000;             % Number of time steps.
dt = mat/T;           % Discrete Time increment dt

%===========================================
% Initialize the terminal stock price matrices 
% for the Euler and Milstein discretization schemes.
%===========================================
Sm = zeros(N,T);          % Milstein
%SE = zeros(N,T);          % Euler stock price
%Se = zeros(N,T);          % Euler log stock price
Sm(:,1) = S0;
%SE(:,1) = S0;
%Se(:,1) = S0;

% Simulate the stock price and the log-stock price under the Euler scheme
% Simulate the stock price under the Milstein schemes
for n=1:N
	for t=2:T
		Z = randn(1);
		W = sqrt(dt)*Z;
		% Euler discretization of stock price
		%SE(n,t) = SE(n,t-1) + (r-q)*SE(n,t-1)*dt + v*SE(n,t-1)*W;
		% Euler discretization of log-stock price
		%Se(n,t) = Se(n,t-1)*exp((r-q-v^2/2)*dt + v*W);
		% Milstein discretization
		Sm(n,t) = Sm(n,t-1) + (r-q)*Sm(n,t-1)*dt + v*Sm(n,t-1)*W...
  			     + 0.5*v^2*Sm(n,t-1)*(Z^2-1)*dt;
        end
    	if strcmp(Averaging,'A')
            A(n) = mean(Sm(n,:));
        elseif strcmp(Averaging,'G')
            A(n) = exp(mean(log(Sm(n,:))));
        end
end

% Terminal stock prices.
%STE = SE(:,end);        % Euler stock price
%STe = Se(:,end);        % Euler log stock price
STm = Sm(:,end);        % Milstein stock price

% Define the payoff
%if strcmp(PutCall,'C')
%	AsianPayoff = max(A - K, 0);
%elseif strcmp(PutCall,'P')
%	AsianPayoff = max(K - A, 0);
%end

% Obtain the simulated option prices.
if strcmp(PutCall,'C')
    % Define the payoff of Vanilla Call.
	BSPrice  = C(S0,K,r,q,v,mat)
   % Define the payoff of Vanilla Call with Milstein scheme.
	MilsteinBS = exp(-r*mat)*mean(max(STm-K,0))
   % Define the payoff of Asian Call.
	AsianPayoff = max(A - K, 0);
else
    % Define the payoff of Vanilla Put.
	BSPrice = P(S0,K,r,q,v,mat)
    % Define the payoff of Vanilla Put with Milstein scheme.
	MilsteinBS = exp(-r*mat)*mean(max(K-STm,0))
    % Define the payoff of Asian Put .
    AsianPayoff = max(K - A, 0);
end
% Calculate the price of the Asian option.
AsianPrice = exp(-r*mat)*mean(AsianPayoff)
Numberofsims = N
NumberofTsteps = T

% Calculate the errors between MilsteinBs Price and Vanilla BsPrice.
%MilsteinError  = abs(MilsteinBS - BSPrice)/BSPrice * 100

