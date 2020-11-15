function example2()

%% Basic parameters
mtrue = [-2.3, 1.2, 2.0, 0.4, -0.8]';

G = rand(7, 5).*3;

% G = [    1,  3.5,    0, -1.5,  1.2;
%          0,  1.5,    2,    2,  0.5;
%       -1.5,  2.3, -0.8,  3.2,  1.4;
%        2.2,  1.4, -1.8,  5.2, -0.4;
%        1.5,  2.8,  0.3, -2.5,  2.9;
%        ];

d = G*mtrue;

%% Conjugate Gradient Method, ART, Kaczmarz and SIRT Method

epsr2norm = 1.0e-6;
epsupdate = 1.0e-5;

% [mkmz, rkmz] = Kaczmarz(G, d, epsr2norm, epsupdate)

% [mart, rart] = ART(G, d, epsr2norm, epsupdate)

% [msrt, rsrt] = SIRT(G, d, epsr2norm, epsupdate)

% [mcgd, rcgd] = ConjugateGradient(G, d, epsr2norm, epsupdate)

end
