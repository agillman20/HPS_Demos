function elliptic_nobody
% this code test the use of the leaf computation with no corner points.
%%%% I have put in sparse matrix operations on the leaf box.  Accuracy
% % % changed to be more accurate.... WHY?? RUN more test


len2    = 1;
Nlevel = 4; % must be divisible by 2;
Npan   = 2^(Nlevel);

% size of a leaf in 1 direction
a       = len2/(2*Npan);
Ncheb = 16;


%%% Pick the equation to be solved
%params  = [1,NaN];               % Laplace problem
%  params  = [2,0.1];               % Modified Helmholtz (i.e. Yukawa)
% params  = [3,2*pi*0.8/(2*a)];    % Helmholtz, kh=params(2)
params  = [4,NaN];               % Variable coefficient Yukawa.
% params  = [5,2*pi*0.8/(2*a)];    % Helmholtz with a "bump" scattering potential.


[NODES] = process_square([0,2*a*Npan,0,2*a*Npan],a,Ncheb,params);


xxs   = [-2;0];
xx_ext = NODES{13,1};
[u,dud1,dud2] = make_ref_sol(xx_ext,params,xxs);
D2N = NODES{23,1};
n = size(D2N,1)/4;
if (mod(n,(Ncheb-2))>0)
    keyboard
end

% keyboard
%%% Test the operators against some known solutions.
num_pts = size(xx_ext,2)/4;
ind_1 = 1:num_pts;
ind_2 = num_pts+1:2*num_pts;
ind_3 = 2*num_pts+1:3*num_pts;
ind_4 = 3*num_pts+1:4*num_pts;
dun = [dud2(ind_1); ...
    dud1(ind_2); ...
    dud2(ind_3); ...
    dud1(ind_4)];
fprintf(1,'max|dun - D2N*u| = %17.8e\n',max(abs(dun - D2N*u)));

[u_ext,dud1_ext,dud2_ext] = create_bdy_data(xx_ext,params,xxs);

leaf_list = makeleaf_list(NODES);

% make a reference solution for a point in the interior
ptlist = [0.25;0];%[0.5,0.2;0.2,0.5];
[uex,~,~] = LOCAL_exact_solution(ptlist,params,xxs);

[u_leaves] = downward_sweep(NODES,u_ext);

[uinterp] = interp_solution(NODES,leaf_list,u_leaves,ptlist,Ncheb);

norm(uinterp- uex)

norm(uinterp- uex)/norm(uex)


keyboard
    

return

function leaf_list = makeleaf_list(NODES)

% this function makes a list of all the leaf boxes

nboxes = size(NODES,2);
leaf_list =[];
for ibox = 1:nboxes
    if NODES{5,ibox} == 0
        leaf_list = [leaf_list, ibox];
    end
    
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [u,dud1,dud2] = create_bdy_data(xx,params,xxs)

if (params(1) == 1)
    
    beta = -1/(2*pi);
    zz1  = xx(1,:)' - xxs(1);
    zz2  = xx(2,:)' - xxs(2);
    u    = 0.5*beta*  log(zz1.*zz1 + zz2.*zz2);
    dud1 =     beta*zz1./(zz1.*zz1 + zz2.*zz2);
    dud2 =     beta*zz2./(zz1.*zz1 + zz2.*zz2);
    
elseif (params(1) == 2)
    
    zz1  = xx(1,:)' - xxs(1);
    zz2  = xx(2,:)' - xxs(2);
    rr   = sqrt(zz1.*zz1 + zz2.*zz2);
    kh   = params(2);
    u    = besselk(0,kh*rr);
    dud1 = -kh*besselk(-1,kh*rr).*zz1./rr;
    dud2 = -kh*besselk(-1,kh*rr).*zz2./rr;
    
elseif (params(1) == 3)
    
    zz1  = xx(1,:)' - xxs(1);
    zz2  = xx(2,:)' - xxs(2);
    rr   = sqrt(zz1.*zz1 + zz2.*zz2);
    kh   = params(2);
    u    = bessely(0,kh*rr);
    dud1 = kh*bessely(-1,kh*rr).*zz1./rr;
    dud2 = kh*bessely(-1,kh*rr).*zz2./rr;
    
elseif (params(1) == 4)
    
    zz1  = xx(1,:)';
    zz2  = xx(2,:)';
    u    = cos(zz1.*zz2);
    dud1 = -zz2.*sin(zz1.*zz2);
    dud2 = -zz1.*sin(zz1.*zz2);
    
else
    
    fprintf(1,'WARNING: No exact solution is coded for this option.\n')
    n    = size(xx,2);
    u    = ones(n,1);
    dud1 = NaN*ones(n,1);
    dud2 = NaN*ones(n,1);
    
end

return

function [u,dud1,dud2] = LOCAL_exact_solution(xx,params,xxs)

if (params(1) == 1)
    
    beta = -1/(2*pi);
    zz1  = xx(1,:)' - xxs(1);
    zz2  = xx(2,:)' - xxs(2);
    u    = 0.5*beta*  log(zz1.*zz1 + zz2.*zz2);
    dud1 =     beta*zz1./(zz1.*zz1 + zz2.*zz2);
    dud2 =     beta*zz2./(zz1.*zz1 + zz2.*zz2);
    
elseif (params(1) == 2)
    
    zz1  = xx(1,:)' - xxs(1);
    zz2  = xx(2,:)' - xxs(2);
    rr   = sqrt(zz1.*zz1 + zz2.*zz2);
    kh   = params(2);
    u    = besselk(0,kh*rr);
    dud1 = -kh*besselk(-1,kh*rr).*zz1./rr;
    dud2 = -kh*besselk(-1,kh*rr).*zz2./rr;
    
elseif (params(1) == 3)
    
    zz1  = xx(1,:)' - xxs(1);
    zz2  = xx(2,:)' - xxs(2);
    rr   = sqrt(zz1.*zz1 + zz2.*zz2);
    kh   = params(2);
    u    = bessely(0,kh*rr);
    dud1 = kh*bessely(-1,kh*rr).*zz1./rr;
    dud2 = kh*bessely(-1,kh*rr).*zz2./rr;
    
elseif (params(1) == 4)
    
    zz1  = xx(1,:)';
    zz2  = xx(2,:)';
    u    = cos(zz1.*zz2);
    dud1 = -zz2.*sin(zz1.*zz2);
    dud2 = -zz1.*sin(zz1.*zz2);
    
else
    
    fprintf(1,'WARNING: No exact solution is coded for this option.\n')
    n    = size(xx,2);
    u    = ones(n,1);
    dud1 = NaN*ones(n,1);
    dud2 = NaN*ones(n,1);
    
end

return




function [NODES] = process_square(box_geom, a,Ncheb,params)


Npan1 = round((box_geom(2) - box_geom(1))/(2*a));
Npan2 = round((box_geom(4) - box_geom(3))/(2*a));

if ( (abs((box_geom(2) - box_geom(1)) - 2*a*Npan1) > 1e-13) || ...
        (abs((box_geom(4) - box_geom(3)) - 2*a*Npan2) > 1e-13) )
    fprintf(1,'ERROR: The box cannot be tesselated into squares of requested size.\n')
    keyboard
end

[LEAF_info,hmin] = construct_leafstuff(Ncheb,a);


% make tree structure.
NODES = make_tree_structure(box_geom,2*a);

% make solver
nboxes = size(NODES,2);
for ibox = nboxes:-1:1
   if NODES{5,ibox}==0 % leaf box
       NODES = process_leaf_box(LEAF_info,NODES,ibox,params,Ncheb);
       
   else
       % decide if it is a vertical or horizontal merge
         if (mod(NODES{2,ibox},2)==0) %horizontal merge
            NODES = LOCAL_mergetwo_hori(NODES,ibox);
        else
            NODES = LOCAL_mergetwo_vert(NODES,ibox);
        end
   end
end


return 
function [u,dud1,dud2] = make_ref_sol(xx,params,xxs)

if (params(1) == 1)
    
    beta = -1/(2*pi);
    zz1  = xx(1,:)' - xxs(1);
    zz2  = xx(2,:)' - xxs(2);
    u    = 0.5*beta*  log(zz1.*zz1 + zz2.*zz2);
    dud1 =     beta*zz1./(zz1.*zz1 + zz2.*zz2);
    dud2 =     beta*zz2./(zz1.*zz1 + zz2.*zz2);
    
elseif (params(1) == 2)
    
    zz1  = xx(1,:)' - xxs(1);
    zz2  = xx(2,:)' - xxs(2);
    rr   = sqrt(zz1.*zz1 + zz2.*zz2);
    kh   = params(2);
    u    = besselk(0,kh*rr);
    dud1 = -kh*besselk(-1,kh*rr).*zz1./rr;
    dud2 = -kh*besselk(-1,kh*rr).*zz2./rr;
    
elseif (params(1) == 3)
    
    zz1  = xx(1,:)' - xxs(1);
    zz2  = xx(2,:)' - xxs(2);
    rr   = sqrt(zz1.*zz1 + zz2.*zz2);
    kh   = params(2);
    u    = bessely(0,kh*rr);
    dud1 = kh*bessely(-1,kh*rr).*zz1./rr;
    dud2 = kh*bessely(-1,kh*rr).*zz2./rr;
    
elseif (params(1) == 4)
    
    zz1  = xx(1,:)';
    zz2  = xx(2,:)';
    u    = cos(zz1.*zz2);
    dud1 = -zz2.*sin(zz1.*zz2);
    dud2 = -zz1.*sin(zz1.*zz2);
    
else
    
    fprintf(1,'WARNING: No exact solution is coded for this option.\n')
    n    = size(xx,2);
    u    = ones(n,1);
    dud1 = NaN*ones(n,1);
    dud2 = NaN*ones(n,1);
    
end

return

function NODES = make_tree_structure(box_geom,lenleaf)

% BOXES is a temporary work array for setting up the tree-structure
BOXES       = zeros(18,100);
lenBOXES    = 100;
BOXES( 1,1) = NaN;
BOXES( 2,1) = 0;
BOXES( 3,1) = NaN;
BOXES(10,1) = NaN;
BOXES(11,1) = box_geom(1);
BOXES(12,1) = box_geom(2);
BOXES(13,1) = box_geom(3);
BOXES(14,1) = box_geom(4);
BOXES(15,1) = -1;
BOXES(16,1) = -1;
BOXES(17,1) = -1;
BOXES(18,1) = -1;
% Create the tree structure by splitting any
% box that holds more than nmax nodes.
%
% We create the boxes one level at a time.
% The following loop is over LEVELS.
ibox_last  = 0;
ibox_new   = 1;
ilevel     = 0;
while (ibox_new > ibox_last)
    
    ibox_first = ibox_last+1;
    ibox_last  = ibox_new;
    ilevel     = ilevel+1;
    
    % Loop over all boxes on the level that was last created.
    % All newly created boxes are temporarily stored in the array TMPBOXES.
    for ibox = ibox_first:ibox_last
        
        % If ibox is larger than lenleaf x lenleaf, it will be partitioned.
        x1min  = BOXES(11,ibox);
        x1max  = BOXES(12,ibox);
        x2min  = BOXES(13,ibox);
        x2max  = BOXES(14,ibox);
        m1     = round((x1max - x1min)/lenleaf); % The horizontal size of the box (counted in nr of leaves).
        m2     = round((x2max - x2min)/lenleaf); % The vertical   size of the box (counted in nr of leaves).
        if (1 < max([m1,m2]))
            
%             indloc = INDS{ibox};
            if (m2 <= m1) % Make a vertical cut.
                m1_west       = round(0.5*m1 + 0.1*(1-2*rand(1))); % How many leaves to put on the left.
                x1half        = x1min + lenleaf*m1_west;
                box_geom_son1 = [x1min,x1half,x2min,x2max];
                box_geom_son2 = [x1half,x1max,x2min,x2max];
            else          % Make a horizontal cut
                m2_south      = round(0.5*m2 + 0.1*(1-2*rand(1))); % How many leaves to put in the south.
                x2half        = x2min + lenleaf*m2_south;
                box_geom_son1 = [x1min,x1max,x2min,x2half];
                box_geom_son2 = [x1min,x1max,x2half,x2max];
            end
            
            % If there is not enough space to save the 2 children in
            % the array BOXES, then double the size of BOXES.
            if ((ibox_new + 2) > lenBOXES)
                BOXES  = [BOXES,zeros(size(BOXES,1),4+size(BOXES,2))];
                lenBOXES = size(BOXES,2);
            end

%             make first child
                ibox_new = ibox_new + 1;
                BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,NaN,...
                    NaN,NaN,NaN,NaN,box_geom_son1,...
                    -1,-1,-1,-1]';
                BOXES(15,ibox) = ibox_new;
%   make second child
                ibox_new = ibox_new + 1;
                BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,NaN,...
                    NaN,NaN,NaN,NaN,box_geom_son2,...
                    -1,-1,-1,-1]';
                BOXES(16,ibox) = ibox_new;
            
        end
        
    end
    
end

nlevels = ilevel - 1;


% Let nboxes denote the number of boxes created,
% create the object NODES which will hold all
% relevant information, and transfer the
% information in BOXES to NODES.
% We also delete the object BOXES.
nboxes = ibox_new;
NODES  = cell(30,nboxes);
for ibox = 1:nboxes
    NODES{01,ibox} = BOXES(11:14,ibox);
    NODES{02,ibox} = BOXES(02,ibox);
    NODES{03,ibox} = BOXES(03,ibox);
    NODES{04,ibox} = [];
    for j = 15:16
        if (BOXES(j,ibox) > 0)
            NODES{04,ibox} = [NODES{04,ibox},BOXES(j,ibox)];
        end
    end
    NODES{05,ibox} = length(NODES{04,ibox});
end

return




function [u] = downward_sweep(NODES,v)

% v = solution on Gaussian nodes on top
%     level.
nboxes = size(NODES,2);

u = cell(2,nboxes);

u{1,1} = v;

%  u(1,:) =  the solution
%      at the gaussian points on the boundary
%      of the box.

% n = Ncheb-2;

for ibox = 1:size(NODES,2)
    if (NODES{5,ibox} > 0)
        %%% ibox is a parent. In this case, map solution on the exterior Gauss
        %%% nodes to the solution on the Gauss nodes on the interface.
%         indext    = NODES{13,ibox};
        indint    = NODES{14,ibox};
        % update children
        ison1 = NODES{4,ibox}(1);
        ison2 = NODES{4,ibox}(2);
        uloc = u{1,ibox};
        vloc = NODES{24,ibox}*u{1,ibox};  % find the solution on the interior for child1
     
    %  propogate solution to children boxes
%         a1    = 0.5*(NODES{01,ison1}(2) - NODES{01,ison1}(1));
        if (mod(NODES{2,ibox},2)==0) % horizontal merge
%             ((NODES{01,ison1}(1) + a1) < NODES{01,ison2}(1)) % horizontal merge
            nloc = size(indint,2);
            nedge = nloc/2;
            u{1,ison1} = [uloc(1:nedge);vloc;uloc(3*nedge+nloc+1:end)];
            u{1,ison2} = [uloc(nedge+1:3*nedge+nloc);vloc(end:-1:1)];
        else
            nedge = size(indint,2);
            u{1,ison1} = [uloc(1:2*nedge);vloc;uloc(5*nedge+1:end)];
            u{1,ison2} = [vloc(end:-1:1);uloc(2*nedge+1:5*nedge)];
        end
    else
        %%% ibox is a leaf. In this case, map solution on the exterior Gauss
        %%% nodes to the solution on the Chebyshev nodes associated with the leaf.
        
        u{2,ibox} = NODES{24,ibox}*u{1,ibox};
        
    end
end

return




function [uinterp] = interp_solution(NODES,leaf_list,u,ptlist,Ncheb)

npts = size(ptlist,2);

nleaves = length(leaf_list);

uinterp = zeros(npts,1);

for j = 1:nleaves
    jbox = leaf_list(j);
    box_geom = NODES{1,jbox};
    xc1 = 0.5*(box_geom(1)+box_geom(2));
    xc2 = 0.5*(box_geom(3)+box_geom(4));
    a = (box_geom(2)-box_geom(1))/2;
    
    %    find the points that live in the box
    
    ind = find((abs(ptlist(1,:)-xc1)<=a).*(abs(ptlist(2,:)-xc2)<=a));
    
    if ~isempty(ind)
        xx = NODES{10,jbox};
        hmin = xx(2,2)-xx(2,1);
        a = 0.5*(box_geom(2)-box_geom(1));
        Jint       = find( max(abs([xx(1,:) - xc1; xx(2,:) - xc2])) < (a - 0.5*hmin));
        uloc = u{2,jbox};
%         Nchebf = sqrt(size(xx,2));
       L = interp2D(xx(:,Jint),ptlist(:,ind),Ncheb-2);
        uinterp(ind,:) = L*uloc(Jint,:);
    end
    
end

return



function L = interp2D(xxin,xout,Ncheb)

% create the vectors for the points in each 
% dimension.


%  First create lagrange interpolants for x direction
 x = xout(1,:); pointx = xxin(1,1:Ncheb-2:end);
 Lx = LOCAL_interpolation_matrix(x,pointx);
% create lagrange interpolant in y direction
 y = xout(2,:); pointy = xxin(2,1:Ncheb-2);
 Ly = LOCAL_interpolation_matrix(y,pointy);
 
ncheb = Ncheb-2;
n = size(xout,2);
nin = size(xxin,2);
L = ones(n,nin);

for i = 1:ncheb
        L(:,(i-1)*ncheb+(1:ncheb)) = (Lx(:,i)*ones(1,ncheb)).*Ly;
end

return

function [LEAF_info,hmin] = construct_leafstuff(Ncheb,a)
%%% Construct the Chebyshev structure for a leaf.
[L,D1,D2,xvec] = LOCAL_get_L_and_D(Ncheb,a);
hmin           = xvec(2) - xvec(1);


Jint = zeros(1,(Ncheb-2)^2);
for j = 2:Ncheb-1
    ind = (j-2)*(Ncheb-2)+1:(j-1)*(Ncheb-2);
    Jint(ind) = (j-1)*Ncheb+2:(j-1)*Ncheb+Ncheb-1;
end


Jext = [     1:Ncheb:Ncheb*(Ncheb-1)...
            Ncheb*(Ncheb-1)+1:Ncheb*Ncheb ...
            Ncheb*(Ncheb-1):-Ncheb:2*Ncheb ...
            Ncheb:-1:2];


indextc= [1 Ncheb 2*Ncheb-1 3*Ncheb-2];

[ZZ1,ZZ2] = meshgrid(xvec);
zz        = [reshape(ZZ1,1,numel(ZZ1));...
    reshape(ZZ2,1,numel(ZZ2))];


indw =Ncheb-1:-1:2;
inds = Ncheb+1:Ncheb:Ncheb*(Ncheb-1);
indn = Ncheb*(Ncheb-1):-Ncheb:2*Ncheb;
inde = Ncheb*(Ncheb-1)+2:Ncheb*Ncheb-1;


LEAF_info = cell(12,1);
LEAF_info{1,1} = L;
LEAF_info{2,1} = D1;
LEAF_info{3,1} = D2;
LEAF_info{4,1} = xvec;
LEAF_info{5,1} = Jint;
LEAF_info{6,1} = Jext;
LEAF_info{7,1} = indextc;
LEAF_info{8,1} = zz;
LEAF_info{9,1} = indw;
LEAF_info{10,1} = inds;
LEAF_info{11,1} = inde;
LEAF_info{12,1} = indn;


return



function L = LOCAL_interpolation_matrix(x,pointx)

n  = size(pointx,2);
L2 = ones(n,size(x,2));

for i=1:n
  for j=1:n
    if (i~=j)
      L2(i,:)=L2(i,:).*(x-pointx(j))/(pointx(i)-pointx(j));
    end
  end
end
L = L2.';  

return

function [NODES] = process_leaf_box(LEAF_info,NODES,ibox,params,Ncheb)
% DtN code
% this leaf procedure removes the corner points. 


% unpack all leaf information
L = LEAF_info{1,1};
D1 =LEAF_info{2,1};
D2 =LEAF_info{3,1};
Jint = LEAF_info{5,1};
Jext= LEAF_info{6,1} ;
indextc = LEAF_info{7,1};
zz = LEAF_info{8,1};
indw= LEAF_info{9,1};
inds = LEAF_info{10,1};
inde = LEAF_info{11,1};
indn = LEAF_info{12,1};

box_geom = NODES{1,ibox};
xc1    = 0.5*(box_geom(1) + box_geom(2));
xc2    = 0.5*(box_geom(3) + box_geom(4));
a      = 0.5*(box_geom(2) - box_geom(1));
len = a;

npts = Ncheb^2;
xx = [xc1;xc2]*ones(1,npts)+len*zz/a;
xx2 = xx;
% list corner points
indcorner = [1,Ncheb, Ncheb^2-Ncheb+1,Ncheb^2];
xx(:,indcorner) = [];

% NODES{10,ibox} = leaf points
NODES{10,ibox} = xx;

zz1 = xx2(1,:)';
zz2 = xx2(2,:)';


Jext2 = Jext;
% Jext2 is all the points on boundary of the original box.
% Jext does not have the corner points anymore.


F = sparse(length(Jext2),npts);
F(:,Jext2) = eye(length(Jext2));
F(:,indcorner) = [];
F(indextc,:) = [];

Jext(indextc) =[];

%%% Construct the full elliptic operator L.
if (params(1) == 1)
    A = -L;
elseif (params(1) == 2)
    kh = params(2);
    A  = -L + kh*kh*eye(size(L,1));
elseif (params(1) == 3)
    kh = params(2);
    A  = -L - kh*kh*eye(size(L,1));
elseif (params(1) == 4)
    b  = zz1.*zz1 + zz2.*zz2;
    A  = -L - diag(b);
elseif (params(1) == 5)
    kh   = params(2);
    b    = 1 - 0.99*exp(-2*(zz1-0.5).^2-2*(zz2-0.5).^2);
    A    = -L - diag(kh*kh*b);
elseif (params(1) == 6)
    kh   = params(2);
    b    = 1 - 0.99*exp(-10*(zz1-0.6).^2-10*(zz2-0.6).^2);
    A    = -L - diag(kh*kh*b);
elseif (params(1) == 7)
    A = -L - 100*D2;
elseif (params(1) == 8)
    kh  = params(2);
    b   = 1-(sin(4*pi*zz1).*sin(4*pi*zz2)).^2;
    A   = -L - kh*kh*diag(b);
elseif (params(1) == 9)
    kh  = params(2);
    b   = 1-exp(-20*(zz1-2).^2 - 20*(zz2-0.5).^2);
    A   = -L - kh*kh*diag(b);
end

A = sparse(A);
Atmp = A(Jint,:);
Atmp(:,indcorner) = [];
D1b = D1;
D1b(:,indcorner) = [];
D2b = D2;
D2b(:,indcorner) = [];

FA = [F;Atmp];
FA_inv = inv(FA);

% create solution operator
X = FA_inv*[eye(4*Ncheb-8); zeros(numel(Jint),4*Ncheb-8)]; % impose f on bdry (but only at Ncheb-1 pts on each edge), PDE inside

% now create differential derivative matrix.

D1bF = D1b*X;
D2bF = D2b*X;

nedge = Ncheb-2;
nbdy = 4*(Ncheb-2);

D = zeros(nbdy);
D(1:nedge,:) = D2bF(inds,:);
D(nedge+1:2*nedge,:) = D1bF(inde,:);
D(2*nedge+1:3*nedge,:) = D2bF(indn,:);
D(3*nedge+1:4*nedge,:) = D1bF(indw,:);


NODES{23,ibox} = D;
NODES{24,ibox} = X;

NODES{13,ibox} = xx2(:,Jext);
NODES{14,ibox} = xx2(:,Jext2);
return


function NODES = LOCAL_mergetwo_hori(NODES,ibox)

%%% Extract the relevant information from "NODES":
ison1 = NODES{4,ibox}(1);
ison2 = NODES{4,ibox}(2);
Tw    = NODES{23,ison1};
Te    = NODES{23,ison2};
%  these are point locations on the boundary of the box
yyw = NODES{13,ison1};
yye = NODES{13,ison2};
nedge = length(yyw)/6;
indw  = 1:length(yyw);
inde  = 1:length(yye);



xxcw  = [mean(NODES{01,ison1}([1,2]));...
         mean(NODES{01,ison1}([3,4]))];
xxce  = [mean(NODES{01,ison2}([1,2]));...
         mean(NODES{01,ison2}([3,4]))];
if (xxce(1) < xxcw(1))
  fprintf(1,'ERROR: The merge assumes "ison1" is the south box.\n')
  keyboard
end


%%% Set up the three index vectors.
J1w = [1:nedge 3*nedge+1:6*nedge];
J3w = nedge + 1:3*nedge;
J2e = [1:4*nedge];
J3e = 6*nedge:-1:4*nedge+1;

%%% Construct the solution operator.
U = (Tw(J3w,J3w) - Te(J3e,J3e))\[-Tw(J3w,J1w), Te(J3e,J2e)];

%%% Construct the new NfD operator;
T = [Tw(J1w,J1w), zeros(length(J1w),length(J2e));...
     zeros(length(J2e),length(J1w)), Te(J2e,J2e)] + ...
    [Tw(J1w,J3w); Te(J2e,J3e)] * U;

%%% Assemble the nodes in the external ring, and order them appropriately.
yyext = [yyw(:,indw(J1w)),yye(:,inde(J2e))];


indtmp = [1:nedge 4*nedge+1:5*nedge 5*nedge+1:8*nedge nedge+1:4*nedge];

%%% Store away various objects 
NODES{13,ibox} = yyext(:,indtmp);    % Index vector for the external (gauss) nodes.
NODES{14,ibox} = yyw(:,J3w);         % Index vector for the external (gauss) nodes.
NODES{23,ibox} = T(indtmp,indtmp);  % NfD operator   [gauss-ext] <- [gauss-ext]
NODES{24,ibox} = U(:,indtmp);       % Solve operator [gauss-int] <- [gauss-ext]


return

function NODES = LOCAL_mergetwo_vert(NODES,ibox)

%%% Extract the relevant information from "NODES":
ison1 = NODES{4,ibox}(1);
ison2 = NODES{4,ibox}(2);
Ts    = NODES{23,ison1};
Tn    = NODES{23,ison2};
yys = NODES{13,ison1};
yyn = NODES{13,ison2};



%%% Extract geometric information about the two boxes.
xxcs  = [mean(NODES{01,ison1}([1,2]));...
    mean(NODES{01,ison1}([3,4]))];
xxcn  = [mean(NODES{01,ison2}([1,2]));...
    mean(NODES{01,ison2}([3,4]))];
if (xxcn(2) < xxcs(2))
    fprintf(1,'ERROR: The merge assumes "ison1" is the south box.\n')
    keyboard
end

nedge = length(yys)/4;

%%% Set up the three index vectors.
J1s = [1:2*nedge 3*nedge+1:4*nedge];
J3s = 2*nedge+1:3*nedge;
J2n = [nedge+1:4*nedge];
J3n = nedge:-1:1;


% %%% Construct the solution operator.
U = (Ts(J3s,J3s) - Tn(J3n,J3n))\[-Ts(J3s,J1s), Tn(J3n,J2n)];

%%% Construct the new NfD operator;
T = [Ts(J1s,J1s), zeros(length(J1s),length(J2n));...
    zeros(length(J2n),length(J1s)), Tn(J2n,J2n)] + ...
    [Ts(J1s,J3s); Tn(J2n,J3n)] * U;

%%% Assemble the nodes in the external ring, and order them appropriately.
yyext = [yys(:,J1s),yyn(:,J2n)];
indtmp = [1:2*nedge 3*nedge+1:6*nedge 2*nedge+1:3*nedge];

% % % %%% Store away various objects
NODES{13,ibox} = yyext(:,indtmp);    % Index vector for the external (gauss) nodes.
NODES{14,ibox} = yys(:,J3s);         % Index vector for the internal (gauss) nodes.
NODES{23,ibox} = T(indtmp,indtmp);  % NfD operator   [gauss-ext] <- [gauss-ext]
NODES{24,ibox} = U(:,indtmp);       % Solve operator [gauss-int] <- [gauss-ext]


return
function [L,D1,D2,xvec] = LOCAL_get_L_and_D(N_side,a)

[D,xvec] = LOCAL_cheb(N_side-1);
xvec     = a*xvec(end:(-1):1);
D        = (1/a)*D;
I        = eye(N_side);
D1       = -kron(D,I);
D2       = -kron(I,D);
Dsq      = D^2;
L        = kron(I,Dsq) + kron(Dsq,I);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Chebyshev nodes on [-1,1].
% It also computes a differentiation operator.
% page 54 of the text
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

