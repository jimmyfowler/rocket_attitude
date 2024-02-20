clear

I = 5000; % Moment of Inertia [m^4]

A = [0, 1; 0, 0];

B = [0; 1/I];


Kp = 1;
Ki = 0;
Kd = 1;

s = tf('s');
C = Kp + Ki/s + Kd*s

% C = pid(Kp,Ki,Kd) % this does the same thing

% bobabooey




