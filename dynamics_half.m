function dx = dynamics_half(t,x,n,Design_Parametersn)

q = x(1:n+1);
dq = x(n+2:end);
%#codegen
tx = transpose(x);
M1 = D(Design_Parametersn,tx);
C1 = C(Design_Parametersn,tx);
G1 = G(Design_Parametersn,tx);
u = zeros(n+1,1);
B1 = 0.01*eye(n,n);%friction
ddq = M1\(-C1*dq-G1+u);%-B1*dq
dx = [dq;ddq];
end
