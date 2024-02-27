function [flag,x_post,P_post,bvec,J_out,augcost,num_nodes,constraint,pos_risk] =...
    mapRiskAverseCvx(num_constrain,y,H,P,R,J_l,x_prior)
% Solves RAPS using optimization approach.
% Computes MAP state estimate using the selected measurements.

Jpminus = P^-1;         % state prior info. matrix
Jrminus = R^-1;         % measurement info. matrix

[m,n] = size(H);
lowerbound = zeros(m,1); % lower bound on measurement selection vector b
upperbound = ones(m,1);  % upper bound on measurement selection vector b

% see intlinprog function documentation for details
diagR  = diag(R);               % diagonal entries of measurement covariance
diagRs = sqrt(diagR);
sH     = diag(diagRs) \ H;      % scale row i by {\sigma_i)^{-1}: diagRs^-1 * H
G = (sH.*sH)';            % G in the paper
G = G(1:num_constrain,:);         % remove unconstrained rows
g = diag(J_l - Jpminus);   % g in the paper
g = g(1:num_constrain,:);         % remove unconstrained rows
option = optimoptions(@linprog,'display','off'); % for output supression

constraint = G*ones(m,1)+diag(Jpminus(1:num_constrain,1:num_constrain));
E_P = chol(Jpminus);    % cholesky decomp. of state prior info. matrix
E_R = chol(Jrminus);    % cholesky decomp. of measurement info. matrix

x_max = x_prior + 200;
x_min = x_prior - 200;

% ensure feasibility
bvec = ones(m,1);          % Initial b to be 1
num_nodes = 0;
flag = true;
% if ~all(G*bvec >= g)% && flag_rapid == true
%     max_lhs = G * bvec;
%     ind = find(g > max_lhs);
%     g(ind) = 0.75 * max_lhs(ind);  % loosen constraint
%     flag = false;
% end

u_len = 0;
ind = [];
if ~all(G*bvec >= g)
    % Add slack variable for soft constraint
    max_lhs = G * bvec;
    ind = find(g > max_lhs);
    u_len = length(ind);
    % G*b + Jp >= (G b_max + Jp) - s
    % G*b + s >= G b_max
    g(ind) = max_lhs(ind);
    flag = false;
end
logicalIndex = zeros(num_constrain,1);
logicalIndex(ind) = 1;
G_o = G(logicalIndex==0,:);
g_o = g(logicalIndex==0);
G_u = G(logicalIndex==1,:);
g_u = g(logicalIndex==1);  
u_coe = 50*ones(1,u_len);
if u_len == 0
    cvx_begin quiet
    variable x(n)
    variable b(m) binary
    variable v(n,m)
    % Objective function
    objf = quad_form(x-x_prior, Jpminus);
    for i=1:m,
        objf = objf + square(y(i)*b(i)-H(i,:)*v(:,i))/diagR(i);
    end
    minimize(objf)
    subject to
    G*b >= g
    x_min <= x <= x_max
    for j=1:n,
        -v(j,:)'+x_max(j)*b >= zeros(m,1);
        v(j,:)'-ones(m,1)*x(j)-x_max(j)*b >= -x_max(j)*ones(m,1);
        v(j,:)'-x_min(j)*b >= zeros(m,1);
        -v(j,:)'+ones(m,1)*x(j)+x_min(j)*b >= x_min(j)*ones(m,1);
    end
    cvx_end
else
    if ~isempty(g_o)
        cvx_begin quiet
        variable x(n)
        variable b(m) binary
        variable v(n,m)
        variable u(u_len)
        % Objective function
        objf = quad_form(x-x_prior, Jpminus)+u_coe*u;
        for i=1:m,
            objf = objf + square(y(i)*b(i)-H(i,:)*v(:,i))/diagR(i);
        end
        minimize(objf)
        subject to
        G_o*b >= g_o
        G_u*b + u >= g_u
        0 <= u <= g_u
        x_min <= x <= x_max
        for j=1:n,
            -v(j,:)'+x_max(j)*b >= zeros(m,1);
            v(j,:)'-ones(m,1)*x(j)-x_max(j)*b >= -x_max(j)*ones(m,1);
            v(j,:)'-x_min(j)*b >= zeros(m,1);
            -v(j,:)'+ones(m,1)*x(j)+x_min(j)*b >= x_min(j)*ones(m,1);
        end
        cvx_end
    else
        cvx_begin quiet
        variable x(n)
        variable b(m) binary
        variable v(n,m)
        variable u(u_len)
        % Objective function
        objf = quad_form(x-x_prior, Jpminus)+u_coe*u;
        for i=1:m,
            objf = objf + square(y(i)*b(i)-H(i,:)*v(:,i))/diagR(i);
        end
        minimize(objf)
        subject to
        G_u*b + u >= g_u
        0 <= u <= g_u
        x_min <= x <= x_max
        for j=1:n,
            -v(j,:)'+x_max(j)*b >= zeros(m,1);
            v(j,:)'-ones(m,1)*x(j)-x_max(j)*b >= -x_max(j)*ones(m,1);
            v(j,:)'-x_min(j)*b >= zeros(m,1);
            -v(j,:)'+ones(m,1)*x(j)+x_min(j)*b >= x_min(j)*ones(m,1);
        end
        cvx_end
    end
end
% augcost = cost(x,x_prior,Jpminus,b,H,y,Jrminus);
augcost = objf;
pos_risk = compute_pos_risk(x,x_prior,P,b(1:m/2),...
    H(1:m/2,:),y(1:m/2),Jrminus(1:m/2,1:m/2));
J_out   = calcJb(b,H,R,Jpminus); % posterior state information matrix
x_post = x;
bvec = b;

%J_out % should meet spec
P_post = inv(J_out);

end

function [C] = cost(x,prior,Jpminus,b,H,y,Jrminus)
ex = x-prior;
ey = y - H*x;
Pb = diag(b);
if isempty(b)
    C1 = ex' * Jpminus * ex;
    C2 =  ey' *  Jrminus *  ey;
    C  = C1 + C2;
else
    C = ex' * Jpminus * ex + ey' * Pb * Jrminus * Pb * ey;
end
end

function [C] = compute_pos_risk(x,prior,P,b_pos,H_pos,y_pos,Jrminus)
ex = x-prior;
ex(4:9) = [];
ex(end) = [];
P_pos = P;
P_pos(4:9,:) = [];
P_pos(:,4:9) = [];
P_pos(:,end) = [];
P_pos(end,:) = [];
Pi = inv(P_pos);
ey = y_pos - H_pos*x;
Pb = diag(b_pos);
C = ex' * Pi * ex + ey' * Pb * Jrminus * Pb * ey;
end

%--------------------------------------------------------------------------
function [x_post,aug_cost] = MAP(by,y,H,E_R,E_P,x_prior)
% Maximum A Posteriori state estimate and the augmented cost function
% when doing measurement selection
% INPUT: by  - measurement selection vector
%        E_P - sqrt(P^-1)
%        E_R - sqrt(R^-1)

if isempty(by)
    error('Input arg "by" is empty.');
else
    Phiby = diag(by);
    A = [E_R * Phiby * H; E_P];             % eqn. (10) in [1]
    c = [E_R * Phiby * y; E_P * x_prior];   % eqn. (10) in [1]
    x_post = (A'*A)^-1*A'*c; % LeastSquares % posterior state estimate
    aug_cost = norm(A*x_post-c)^2;          % eqn. (11) in [1]
end
end
%--------------------------------------------------------------------------
function J = calcJb(by,H,R,Jpminus)
% Calculates posterior state information matrix when doing measurement selection
% INPUT: by - measurement selection vector
%        Jpminus - State information matrix prior, Jpminus = Pminus^-1

Pby  = diag(by); % Phi(by)
PhiH = Pby * H;
J    = PhiH' * R^-1 * PhiH + Jpminus; % eqn. (13) in [1]
end