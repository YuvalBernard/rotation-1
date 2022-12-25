function [Z,A,domain] = CEST_numerical(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,db,w1,tp)

R1a = 1/T1a; R2a = 1/T2a;
R1b = 1/T1b; R2b = 1/T2b;
ka = kb*M0b/M0a; % Given from mass equality on both sides of exchange chemical equation
M0 = [0 0 0 0 M0a M0b 1]';
dwb = dwa + db;
Z = zeros(size(dwa));

for j = 1:length(dwa)
    Mprime = @(y) [
        -(R2a+ka)*y(1)+kb*y(2)+dwa(j)*y(3);
        ka*y(1)-(R2b+kb)*y(2)+dwb(j)*y(4);
        -dwa(j)*y(1)-(R2a+ka)*y(3)+kb*y(4)+w1*y(5);
        -dwb(j)*y(2)+ka*y(3)-(R2b+kb)*y(4)+w1*y(6);
        -w1*y(3)-(R1a+ka)*y(5)+kb*y(6)+R1a*M0a*y(7);
        -w1*y(4)+ka*y(5)-(R1b+kb)*y(6)+R1b*M0b*y(7);
        0
        ];
%     Mprime = @(y) A*y;
    [~,y] = ode45(@(t,y) Mprime(y),[0 tp],M0);    
    Z(j) = y(end,5)/M0a;
end
%  Assuming the domain is symmetric
domain = round(length(dwa)/2);
A = fliplr(Z(domain:end)) - Z(1:domain);