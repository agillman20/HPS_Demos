function main_1D_HPS
% this is a full 1D HPS code

% number of chebyshev discretization points per leaf
Ncheb = 16;
%box_geom = geometry for which the problem is considered
%box_geom = domain = Omega
box_geom = [0 1];

%nlevel = number of levels in the tree;
nlevel = 5;

% precompute the solver
NODES = create_solver(Ncheb,box_geom,nlevel);

xxext = NODES{13,1};
% create the boundary data
[ubdy,~] = create_bdy_data(xxext);

% build the particular solution operators
[u] = upward_sweep(NODES);

% Create solution that specified by boundary 
% data

[u] = downward_sweep(NODES,u,ubdy');


% Compute the error
leaf_list = makeleaf_list(NODES);

meanrel = 0;
for j = 1:length(leaf_list)
    jbox = leaf_list(j);
    xx = NODES{10,jbox};
%     uex = sin(xx);
alpha = 1/(5*pi);
uex = sin(1./(alpha+xx));

meanrel = meanrel+norm(u{2,jbox}-uex)/norm(uex);
end
meanrel/length(leaf_list)
(Ncheb-1)*2^nlevel

keyboard

return

function NODES = create_solver(Ncheb,box_geom,nlevel)

% create the binary tree
NODES = create_tree(box_geom,nlevel);

nboxes = size(NODES,2);

for ibox = nboxes:-1:1
    if NODES{5,ibox} ==0
        NODES = process_leaf(ibox,NODES,Ncheb);
    else
        NODES = mergetwo(NODES,ibox);
    end
    
    
end



return


function NODES = create_tree(box_geom,nlevel)

nboxes = 0;
for j = 0:nlevel
    nboxes = nboxes+2^j;
end

NODES = cell(50,nboxes);

NODES{1,1} = box_geom;
NODES{2,1} = 0;

icount = 1;
for ibox = 1:nboxes-2^nlevel
    ison1 = icount+1;
    ison2 = icount+2;
    box_geom = NODES{1,ibox};
    mid = 0.5*(box_geom(1)+box_geom(2));
    NODES{1,ison1} = [box_geom(1) mid];
    NODES{1,ison2} = [mid box_geom(2)];
    NODES{4,ibox} = [ison1 ison2];
    NODES{5,ibox} = 2;
    NODES{3,ison1} = ibox;
    NODES{3,ison2} = ibox;
    NODES{2,ison1} = NODES{1,ibox}+1;
    NODES{2,ison2} = NODES{1,ibox}+1;
    NODES{5,ison1} = 0;
    NODES{5,ison2} = 0;
    icount = icount+2;
end

return


function src = create_source(xx)
% xx are the interior points
% src = -sin(xx);

alpha = 1/(5*pi);
src  = -2./(alpha+xx).^3.*cos(1./(alpha+xx));

return

function [ubdy,un] = create_bdy_data(xx)
% xx are the points on the boundary

% % u at boundary
% ubdy = sin(xx);
% % flux of solution at boundary
% un = cos(xx);

alpha = 1/(5*pi);
ubdy = sin(1./(alpha+xx));
un = -1./(alpha+xx).^2.*cos(1./(alpha+xx));


return

function [xvec,D] = get_cheb(N_side,a)

[D,xvec] = LOCAL_cheb(N_side-1);
xvec     = a*xvec(end:(-1):1);
D        = (1/a)*D;


return

function [NODES] = process_leaf(ibox,NODES,Ncheb)

box_geom = NODES{1,ibox};
a = 0.5*(box_geom(2)-box_geom(1));
xc = 0.5*(box_geom(2)+box_geom(1));

[xvec,D] = get_cheb(Ncheb,a);


xx = xc + xvec;
hmin = xvec(2)-xvec(1);

% get ordering for points.
Jint     = find(abs(xvec) < (a - 0.5*hmin));
Jext = LOCAL_complement_vector(Jint,Ncheb);

%create the derivative operator. Here it d^2/dx^2
% A = D^2;

alpha = 1/(5*pi);
d = 1./(alpha+xx).^4;
A = -D^2-diag(d);


F                = inv(A(Jint,Jint));
% Create particular solution operator
S             = zeros(Ncheb,length(Jext));
S(Jext,:) = eye(length(Jext));
S(Jint,:)     =  -F*A(Jint,Jext);

NODES{24,ibox} = S;
T  = -D(Jext,:)*S;
NODES{23,ibox} = T;
% create homogenous solution operator.
Fbig         = zeros(Ncheb,length(Jint));
Fbig(Jint,:) = F;
H = -D(Jext,:)*Fbig;

% particular solution
%  the solution operator
NODES{25,ibox} = Fbig;
%  for creating normal derivative.
NODES{26,ibox} = H;

NODES{10,ibox} = xx;
NODES{13,ibox} = xx(Jext);
NODES{28,ibox} = Jint;


return


function indc = LOCAL_complement_vector(indi,N)
indtmp = 1:N;
indtmp(indi) = 2*N*ones(1,length(indi));
indtmp = sort(indtmp);
indc = indtmp(1:(N - length(indi)));
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Chebyshev nodes on [-1,1].
% It also computes a differentiation operator.
% It is modified from a code by Nick Trefethen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D,x] = LOCAL_cheb(N)

if N==0
    D=0;
    x=1;
    return
end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
D  = D - diag(sum(D'));                 % diagonal entries

return

function NODES = mergetwo(NODES,ibox)

ison1 = NODES{4,ibox}(1);
ison2 = NODES{4,ibox}(2);
Tw    = NODES{23,ison1};
Te    = NODES{23,ison2};
yyw = NODES{13,ison1};
yye = NODES{13,ison2};

NODES{1,ibox} = [NODES{1,ison1}(1) NODES{1,ison2}(2)];

%%% Set up the three index vectors.
J1w = find(yyw == (NODES{01,ison1}(1)));
J3w = find(yyw > (NODES{01,ison1}(1)));
J2e = find(yye == (NODES{01,ison2}(2)));
J3e = find(yye < (NODES{01,ison2}(2) ));

% %%% Construct the solution operator.
U = (Tw(J3w,J3w) - Te(J3e,J3e))\[-Tw(J3w,J1w), Te(J3e,J2e)];

%%% Construct the new NfD operator;
T = [Tw(J1w,J1w), zeros(length(J1w),length(J2e));...
    zeros(length(J2e),length(J1w)), Te(J2e,J2e)] + ...
    [Tw(J1w,J3w); Te(J2e,J3e)] * U;

%%% Assemble the nodes in the external ring, and order them appropriately.

yyext = [yyw(J1w),yye(J2e)];

% % % %%% Store away various objects
NODES{13,ibox} = yyext;    % Index vector for the external (gauss) nodes.
NODES{14,ibox} = yyw(J3w);         % Index vector for the internal (gauss) nodes.
NODES{23,ibox} = T;  % NfD operator   [gauss-ext] <- [gauss-ext]
NODES{24,ibox} = U;       % Solve operator [gauss-int] <- [gauss-ext]
NODES{25,ibox} = inv((Tw(J3w,J3w) - Te(J3e,J3e)));
tmp = [Tw(J1w,J3w); Te(J2e,J3e)]*NODES{25,ibox};
NODES{26,ibox} =tmp;% tmp(indtmp,:);
NODES{28,ibox} = J3e;


return

function[u] = upward_sweep(NODES)
% v = solution on Gaussian nodes on top
%     level.
nboxes = size(NODES,2);

u = cell(4,nboxes);


% v_n = neumann data

%  u(1,:) =  the solution
%      at the gaussian points on the boundary
%      of the box.

% u{2,:} = particular solution on the boundary

% n = Ncheb-2;

% u{3,:} = neumann data from particular solution


% upward pass ==== collect all particular solutions
% and neumann data

for ibox = nboxes:-1:1
    if NODES{5,ibox} == 0
        xxloc = NODES{10,ibox};
        ind = NODES{28,ibox};
        % compute source
        g = create_source(xxloc(ind));
        %compute solution due to source
        u{2,ibox} = NODES{25,ibox}*g;
        % compute neumann data from particular solution
        u{3,ibox} = NODES{26,ibox}*g;
        
    else
        ison1 = NODES{4,ibox}(1);
        ison2 = NODES{4,ibox}(2);
        hA = u{3,ison1};
        hB = u{3,ison2};
        
        h3A = hA(2);%#ok
        h3B = hB(1);%#ok
        
        u{2,ibox} = NODES{25,ibox}*(-h3A+h3B);
        
        he = [hA(1);hB(2)] + NODES{26,ibox}*(-h3A+h3B);
        u{3,ibox} = he;
        
    end
    
end

return

function [u] = downward_sweep(NODES,u,v)

% solution on the boundary
u{1,1} = v;


for ibox = 1:size(NODES,2)
    if (NODES{5,ibox} > 0)
        %%% ibox is a parent. In this case, map solution on the exterior Gauss
        %%% nodes to the solution on the Gauss nodes on the interface.
        % update children
        ison1 = NODES{4,ibox}(1);
        ison2 = NODES{4,ibox}(2);
        uloc = u{1,ibox};
        vloc = NODES{24,ibox}*u{1,ibox}+u{2,ibox};  % find the solution on the interior for child1
        
        %  propogate solution to children boxes
        u{1,ison1} = [uloc(1);vloc];
        u{1,ison2} = [vloc;uloc(2)];
        
    else
        %%% ibox is a leaf. In this case, map solution on the exterior Gauss
        %%% nodes to the solution on the Chebyshev nodes associated with the leaf.
        %                     Homogenous solution   + particular
        u{2,ibox} = NODES{24,ibox}*u{1,ibox}+u{2,ibox};
        
    end
end

return



function leaf_list = makeleaf_list(NODES)

% this function makes a list of all the leaf boxes

nboxes = size(NODES,2);
leaf_list =[];
for ibox = 1:nboxes
    if NODES{5,ibox} == 0;
        leaf_list = [leaf_list, ibox];
    end
    
end
return
