function L = findKL(KL, Input,B,C,D,G)

desired_omega = [Input.omega];
Ln = KL(end-Input.n:end-1); %length of each link
Lcenter = KL(end);
L = 0;

L_total(1) = Lcenter;
L_total= [L_total Ln];
[mass,  COM_x,  COM_y, Ixx, Iyy]  = gen_SASA_model_params_poly(Input.X,Input.Y,Input.t,L_total);

Mn = mass(2:end); %mass of each link
Jn = Iyy(2:end); %inertia of each link
MR_COM_xn = COM_x(2:end); % COM position for each link
MR_COM_yn = COM_y(2:end);

cd SimMatrixFun
for i=1:size(Input.x,2)
    % M1 = D([Input.DParams(1:Input.n), Ln ,Input.DParams((2*Input.n+1):end-1),Lcenter],transpose(Input.x(:,i)));
    % M1 = feval(strcat('D',num2str(Input.n)),[Mn, Ln ,Input.DParams((2*Input.n+1):end-1),Lcenter],transpose(Input.x(:,i)));
    %i
    % [Mn, Ln ,KL(i)*ones(1,(Input.n)),MR_COM_xn,MR_COM_yn,Lcenter],transpose(Input.x(:,i))
    M1 = D([Mn, Ln,Jn ,KL(i)*ones(1,(Input.n)),MR_COM_xn,MR_COM_yn,Lcenter],transpose(Input.x(:,i)));
    K_diag = ones(1,(Input.n));
    for k = 1:(Input.n)
        K_diag(k) = sqrt(205e9*Jn(k)/Mn(k));
    end
    K_diag = KL(i)*K_diag;
    G = diag([K_diag]);
    
    % b = B([Mn, Ln ,KL(i)*ones(1,(Input.n)),MR_COM_xn,MR_COM_yn,Lcenter],transpose(Input.x(:,i)));
    
    A = pinv(M1)*G;
    % A = pinv(M1)*b;
    
    eigA = sort(eig(A));
    % % eigA/norm(eigA)
    % F = eigA(2) - desired_omega^2;
    % F = [(eigA(1:length(desired_omega)) - desired_omega.^2)];
    idx = min(Input.n,length(desired_omega));
        F = [norm((eigA(1:idx) - desired_omega(1:idx).^2))/norm(desired_omega(1:idx).^2)];
    
    L = L+F;
end
L = L/size(Input.x,2)*100;
cd ..
end
