function Laplace_bodyload
% this code test the use of the leaf computation with no corner points.

len2    = 1;
Nlevel = 1; % must be divisible by 2;
Npan   = 2^(Nlevel);

% size of a leaf in 1 direction
a       = len2/(2*Npan);
Ncheb = 16;


%%% Pick the equation to be solved
%params  = [1,NaN];               % Laplace problem
% params  = [2,0.1];               % Modified Helmholtz (i.e. Yukawa)
% params  = [3,2*pi*0.8/(2*a)];    % Helmholtz, kh=params(2)
%params  = [4,NaN];               % Variable coefficient Yukawa.
% params  = [5,2*pi*0.8/(2*a)];    % Helmholtz with a "bump" scattering potential.
params = [10,NaN]; %Poisson with bodyload


[NODES] = Build_sol_wbd([0,2*a*Npan,0,2*a*Npan], a,Ncheb,params);


u = upward_sweep(NODES,Ncheb,params);


% verify the flux
xbdy = NODES{13,1};
xxs   = [-2;0];
[ubdy,ux1,uy1] = create_bdy_data(xbdy,params,xxs);


nedge = size(xbdy,2)/4;

unex = [uy1(1:nedge) ux1(nedge+1:2*nedge) uy1(2*nedge+1:3*nedge) ux1(3*nedge+1:4*nedge)];

[u] = downward_sweep(NODES,u,ubdy,params);


% check on the performance in one leaf box
ibox = 7;

uapp = u{2,ibox};


xx = NODES{10,ibox};
indcorner = [1,Ncheb, Ncheb^2-Ncheb+1,Ncheb^2];
xx2 = xx;
xx2(:,indcorner) = [];


norm(uapp)
[uex,~,~] = LOCAL_exact_solution(xx2,params,xxs);
norm(uex)

norm(uapp-uex)
norm(uapp-uex)/norm(uex)

keyboard



return

function g = create_source(params,xx)

n = size(xx,2);
if params(1)<10
    g = zeros(n,1);
else
    if params(1) == 10
    alpha = 1;
    xxc = 0.51;
    yyc = 0.117;
    zz1 = xx(1,:)';
    zz2 = xx(2,:)';
    u = exp(-alpha*((zz1-xxc).^2+(zz2-yyc).^2));
%     ux = -2*alpha*(zz1-xxc).*u;
%     uy = -2*alpha*(zz2-yyc).*u;
    uxx = (4*alpha^2*(zz1-xxc).^2-2*alpha).*u;
    uyy = (4*alpha^2*(zz2-yyc).^2-2*alpha).*u;
        g  = -(uxx+uyy);
    elseif params(1) == 11
%         eps = 1;
    eps = params(2) ;
    zz1 = xx(1,:)';
    zz2 = xx(2,:)';
ux = -((exp((-1 + zz1)/eps).*(1 - exp((-1 + zz2)/eps)).*cos(pi*(zz1 + zz2)))/eps) - ...
      (1 - exp((-1 + zz1)/eps)).*(1 - exp((-1 + zz2)/eps)).*pi.*sin(pi*(zz1 + zz2));

 
uy = -((exp((-1 + zz2)/eps).*(1 - exp((-1 + zz1)/eps)).*cos(pi*(zz1 + zz2)))/eps) - ...
  (1 - exp((-1 + zz1)/eps)).*(1 - exp((-1 + zz2)/eps)).*pi.*sin(pi*(zz1 + zz2)) ;
     
     
uxx =         (-1 + exp((-1 + zz1)/eps)).*(1 - exp((-1 + zz2)/eps)).*pi^2.*cos(pi*(zz1 + zz2)) - ...
 (exp((-1 + zz1)/eps).*(1 - exp((-1 + zz2)/eps)).*cos(pi*(zz1 + zz2)))/eps^2 + ...
 (2*exp((-1 + zz1)/eps).*(1 - exp((-1 + zz2)/eps)).*pi.*sin(pi*(zz1 + zz2)))/eps;
     
     
uyy =         (-1 + exp((-1 + zz1)/eps)).*(1 - exp((-1 + zz2)/eps)).*pi.^2.*cos(pi*(zz1 + zz2)) - ...
 (exp((-1 + zz2)/eps).*(1 - exp((-1 + zz1)/eps)).*cos(pi*(zz1 + zz2)))/eps^2 + ...
 (2*exp((-1 + zz2)/eps).*(1 - exp((-1 + zz1)/eps)).*pi.*sin(pi*(zz1 + zz2)))/eps;         
        
 
 g = -eps*(uxx+uyy)+2*ux+uy;
    end
end

return



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

elseif (params(1) == 10)
    alpha = 1;
    xxc = 0.51;
    yyc = 0.117;
    zz1 = xx(1,:)';
    zz2 = xx(2,:)';
    u = exp(-alpha*((zz1-xxc).^2+(zz2-yyc).^2));
    ux = -2*alpha*(zz1-xxc).*u;
    uy = -2*alpha*(zz2-yyc).*u;
    dud1 = ux;
    dud2 = uy;
    
elseif (params(1) == 11)
    eps = params(2);
        zz1 = xx(1,:)';
    zz2 = xx(2,:)';
ux = -((exp((-1 + zz1)/eps).*(1 - exp((-1 + zz2)/eps)).*cos(pi*(zz1 + zz2)))/eps) - ...
      (1 - exp((-1 + zz1)/eps)).*(1 - exp((-1 + zz2)/eps)).*pi.*sin(pi*(zz1 + zz2));

 
uy = -((exp((-1 + zz2)/eps).*(1 - exp((-1 + zz1)/eps)).*cos(pi*(zz1 + zz2)))/eps) - ...
  (1 - exp((-1 + zz1)/eps)).*(1 - exp((-1 + zz2)/eps)).*pi.*sin(pi*(zz1 + zz2)) ;

 u = (1-exp(-(1-zz1)/eps)).*(1-exp(-(1-zz2)/eps)).*cos(pi*(zz1+zz2));
 
 dud1 = ux;
 dud2 = uy;
 
    
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

elseif (params(1) == 10)
    alpha = 1;
    xxc = 0.51;
    yyc = 0.117;
    zz1 = xx(1,:)';
    zz2 = xx(2,:)';
    u = exp(-alpha*((zz1-xxc).^2+(zz2-yyc).^2));
    ux = -2*alpha*(zz1-xxc).*u;
    uy = -2*alpha*(zz2-yyc).*u;
    dud1 = ux;
    dud2 = uy;
    
elseif (params(1) == 11)
    eps = params(2);
        zz1 = xx(1,:)';
    zz2 = xx(2,:)';
ux = -((exp((-1 + zz1)/eps).*(1 - exp((-1 + zz2)/eps)).*cos(pi*(zz1 + zz2)))/eps) - ...
      (1 - exp((-1 + zz1)/eps)).*(1 - exp((-1 + zz2)/eps)).*pi.*sin(pi*(zz1 + zz2));

 
uy = -((exp((-1 + zz2)/eps).*(1 - exp((-1 + zz1)/eps)).*cos(pi*(zz1 + zz2)))/eps) - ...
  (1 - exp((-1 + zz1)/eps)).*(1 - exp((-1 + zz2)/eps)).*pi.*sin(pi*(zz1 + zz2)) ;

 u = (1-exp(-(1-zz1)/eps)).*(1-exp(-(1-zz2)/eps)).*cos(pi*(zz1+zz2));
 
 dud1 = ux;
 dud2 = uy;
 
    
else
    
    fprintf(1,'WARNING: No exact solution is coded for this option.\n')
    n    = size(xx,2);
    u    = ones(n,1);
    dud1 = NaN*ones(n,1);
    dud2 = NaN*ones(n,1);
    
end

return


function [u] = upward_sweep(NODES,Ncheb,params)
% build particular solution information.

nboxes = size(NODES,2);

u = cell(3,nboxes);
%u(1,:) = solution at the nodes on the boundary 
% of the box
% u(2,:) = particular solution on the boundary
% u(3,:) = flux of particular solution

% upward sweep to make outgoing impedance data he
for ibox = nboxes:-1:1
    if isempty(NODES{4,ibox})
        %        create source data
        x_int = NODES{10,ibox}(:,NODES{28,ibox});
        src = create_source(params,x_int);
        rhs = [zeros(4*Ncheb-4,1);src];
        % flux of particular solution
        u{3,ibox} = NODES{26,ibox}*rhs;

        % make solution contribution from body load
        u{2,ibox} = NODES{27,ibox}*rhs;

    else
        ikid1 = NODES{4,ibox}(1);
        ikid2 = NODES{4,ibox}(2);
        J1 = NODES{31,ikid1};
        J2 = NODES{31,ikid2};
        J3_A = NODES{32,ikid1};
        J3_B = NODES{32,ikid2};

        W = NODES{26,ibox};
        BB = NODES{34,ibox};

        hA = u{3,ikid1};
        hB = u{3,ikid2};

        h3A = hA(J3_A);
        h3B = hB(J3_B);
        indtmp = NODES{33,ibox};

        % particular solution on the J3 
        u{2,ibox} = W*(-h3A+h3B);
        % flux of particular solution on the boundary 
        he = [hA(J1);hB(J2)] +BB*(-h3A+h3B);
        u{3,ibox} = he(indtmp);

    end

end

return

function [u] = downward_sweep(NODES,u,bc,params)

% downward sweep propogating incoming impedance data
% to leaf boxes.  Then solution
nboxes = size(NODES,2);

u{1,1} = bc; % Dirichlet boundary data
for ibox = 1:nboxes

    if isempty(NODES{4,ibox})
        % leaf box
        u{2,ibox} = NODES{24,ibox}*u{1,ibox}+u{2,ibox};
        xx = NODES{10,ibox};
        xxs= [];
    else
        % for non-leaf boxes make incoming impedance data
        ikid1 = NODES{4,ibox}(1);
        ikid2 = NODES{4,ibox}(2);
        vloc1 = NODES{24,ibox}*u{1,ibox}+u{2,ibox};  % find the solution on the interior
        %                                               edge for ison1

        uloc = u{1,ibox};
        indint    = NODES{14,ibox};

        a1    = 0.5*(NODES{01,ikid1}(2) - NODES{01,ikid1}(1));
        if ((NODES{01,ikid1}(1) + a1) < NODES{01,ikid2}(1)) % horizontal merge
            nloc = size(indint,2);
            nedge = nloc/2;
            u{1,ikid1} = [uloc(1:nedge);vloc1;uloc(3*nedge+nloc+1:end)];
            u{1,ikid2} = [uloc(nedge+1:3*nedge+nloc);vloc1(end:-1:1)];
        else
            nedge = size(indint,2);
            u{1,ikid1} = [uloc(1:2*nedge);vloc1;uloc(5*nedge+1:end)];
            u{1,ikid2} = [vloc1(end:-1:1);uloc(2*nedge+1:5*nedge)];
        end

    end

end

return




function [LEAF_info] = construct_leafstuff(Ncheb,a)
%%% Construct the Chebyshev structure for a leaf.
[L,D1,D2,xvec] = LOCAL_get_L_and_D(Ncheb,a);

% create indexing of local nodes
Jint = zeros(1,(Ncheb-2)^2);
Jint2 = Jint;
for j = 2:Ncheb-1
    ind = (j-2)*(Ncheb-2)+1:(j-1)*(Ncheb-2);
    Jint(ind) = (j-1)*Ncheb+2:(j-1)*Ncheb+Ncheb-1;
    Jint2(ind) = (j-1)*Ncheb:j*Ncheb-3;
end

% exterior nodes with respect to tpg with corners
Jext = [     1:Ncheb:Ncheb*(Ncheb-1)...
    Ncheb*(Ncheb-1)+1:Ncheb*Ncheb ...
    Ncheb*(Ncheb-1):-Ncheb:2*Ncheb ...
    Ncheb:-1:2];

indextc= [1 Ncheb 2*Ncheb-1 3*Ncheb-2];


% create tensor product grid
[ZZ1,ZZ2] = meshgrid(xvec);
zz        = [reshape(ZZ1,1,numel(ZZ1));...
    reshape(ZZ2,1,numel(ZZ2))];

% edge point ordering with respect to
% tensor product grid including corners.

indw =Ncheb-1:-1:2;
inds = Ncheb+1:Ncheb:Ncheb*(Ncheb-1);
indn = Ncheb*(Ncheb-1):-Ncheb:2*Ncheb;
inde = Ncheb*(Ncheb-1)+2:Ncheb*Ncheb-1;

% indices for boundary including corner nodes
%%%% This can be removed but requires additional 
%%%% editing of the leaf comptuations
indsl = 1:Ncheb:Ncheb*(Ncheb-1);
indel = Ncheb*(Ncheb-1)+1:Ncheb*Ncheb-1;
indnl = Ncheb*Ncheb:-Ncheb:2*Ncheb;
indwl = Ncheb:-1:2;


LEAF_info = cell(17,1);
LEAF_info{1,1} = L;
LEAF_info{2,1} = D1;
LEAF_info{3,1} = D2;
LEAF_info{5,1} = Jint;
LEAF_info{6,1} = Jext;
LEAF_info{7,1} = indextc;
LEAF_info{8,1} = zz;
LEAF_info{9,1} = indw;
LEAF_info{10,1} = inds;
LEAF_info{11,1} = inde;
LEAF_info{12,1} = indn;
LEAF_info{13,1} = indwl;
LEAF_info{14,1} = indsl;
LEAF_info{15,1} = indel;
LEAF_info{16,1} = indnl;
LEAF_info{17,1} = Jint2;


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NODES = make_tree_structure(box_geom,lenleaf)

% BOXES is a temporary work array for setting up the tree-structure
BOXES       = zeros(18,100);
lenBOXES    = 100;
BOXES( 1,1) = NaN;
BOXES( 2,1) = 0;
BOXES( 3,1) = NaN;
% BOXES( 6,1) = 1;
% BOXES( 7,1) = ntot;
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
                %                 J_son1        = find(xx(1,indloc) <= (x1half + 0.5*hmin));
                %                 J_son2        = find(xx(1,indloc) >= (x1half - 0.5*hmin));
                box_geom_son1 = [x1min,x1half,x2min,x2max];
                box_geom_son2 = [x1half,x1max,x2min,x2max];
            else          % Make a horizontal cut
                m2_south      = round(0.5*m2 + 0.1*(1-2*rand(1))); % How many leaves to put in the south.
                x2half        = x2min + lenleaf*m2_south;
                %                 J_son1        = find(xx(2,indloc) <= (x2half + 0.5*hmin));
                %                 J_son2        = find(xx(2,indloc) >= (x2half - 0.5*hmin));
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
            %             if ~isempty(J_son1)
            ibox_new = ibox_new + 1;
            BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,NaN,...
                NaN,NaN,NaN,NaN,box_geom_son1,...
                -1,-1,-1,-1]';
            BOXES(15,ibox) = ibox_new;
            ibox_new = ibox_new + 1;
            BOXES(:,ibox_new) = [NaN,ilevel,ibox,NaN,NaN,NaN,...
                NaN,NaN,NaN,NaN,box_geom_son2,...
                -1,-1,-1,-1]';
            BOXES(16,ibox) = ibox_new;
            %                 INDS{ibox_new} = indloc(J_son2);
            %             end

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

function [NODES] = Build_sol_wbd(box_geom, a,Ncheb,params)


Npan1 = round((box_geom(2) - box_geom(1))/(2*a));
Npan2 = round((box_geom(4) - box_geom(3))/(2*a));

if ( (abs((box_geom(2) - box_geom(1)) - 2*a*Npan1) > 1e-13) || ...
        (abs((box_geom(4) - box_geom(3)) - 2*a*Npan2) > 1e-13) )
    fprintf(1,'ERROR: The box cannot be tesselated into squares of requested size.\n')
    keyboard
end

[LEAF_info] = construct_leafstuff(Ncheb,a);


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
indwl = LEAF_info{13,1}; 
indsl = LEAF_info{14,1};
indel = LEAF_info{15,1};
indnl = LEAF_info{16,1};
Jint2 = LEAF_info{17,1};

%local numbering for boundary to eliminate the corners
indbd = [2:Ncheb-1 Ncheb+1:2*Ncheb-2 2*Ncheb:3*Ncheb-3 ...
          3*Ncheb-1:4*Ncheb-4];

box_geom = NODES{1,ibox};
xc1    = 0.5*(box_geom(1) + box_geom(2));
xc2    = 0.5*(box_geom(3) + box_geom(4));
a      = 0.5*(box_geom(2) - box_geom(1));

% % create actual points (needed for the adaptive version of the
% % method.
xx = zz+[xc1;xc2]*ones(1,Ncheb*Ncheb);
% % list corner points

% NODES{10,ibox} = leaf points
NODES{10,ibox} = xx;

zz1 = xx(1,:)';
zz2 = xx(2,:)';

% list corner points
indcorner = [1,Ncheb, Ncheb^2-Ncheb+1,Ncheb^2];

Jext2 = Jext;
Jext2(indextc) =[];
% Jext2 does not have the corner points anymore.


F = sparse(length(Jext),Ncheb^2);
F(:,Jext) = eye(length(Jext));
F(:,indcorner) = [];
F(indextc,:) = [];

Um = eye(Ncheb^2);  Um = Um(Jext,:);


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
elseif (params(1) == 10)
    A = -L;
end

 A = sparse(A);
% 
Vm = [D2(indsl,:); ...
    D1(indel,:); ...
    D2(indnl,:); ...
    D1(indwl,:)];

% % % % create solution operator - homogenous

FA2= [Um; A(Jint,:)];
FA2_inv = inv(FA2);
indtmp = 1:Ncheb^2;
indtmp(indcorner) =[];

X2 = FA2_inv*[eye(4*Ncheb-4); zeros(numel(Jint),4*Ncheb-4)];


% homogenous solution operator
Y2 = X2(indtmp,indbd);
% particular solution operator
W2 = FA2_inv(indtmp,:);


Vm2 = Vm(indbd,indtmp);


% now create differential derivative matrix.


DD2 = Vm2*Y2;


% flux of particular solution 

H2 = Vm2*FA2_inv(indtmp,:);


NODES{10,ibox} = xx;
NODES{13,ibox} = xx(:,Jext2);

NODES{23,ibox} = DD2;
NODES{24,ibox} = Y2;

NODES{27,ibox} = W2;
NODES{26,ibox} = H2;

NODES{28,ibox} = Jint;
NODES{29,ibox} = Jint2;


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

function NODES = LOCAL_mergetwo_hori(NODES,ibox)

%%% Extract the relevant information from "NODES":
ikid1 = NODES{4,ibox}(1);
ikid2 = NODES{4,ibox}(2);
Tw    = NODES{23,ikid1};
Te    = NODES{23,ikid2};
%  these are point locations on the boundary of the box
yyw = NODES{13,ikid1};
yye = NODES{13,ikid2};
nedge = length(yyw)/6;
indw  = 1:length(yyw);
inde  = 1:length(yye);


xxcw  = [mean(NODES{01,ikid1}([1,2]));...
         mean(NODES{01,ikid1}([3,4]))];
xxce  = [mean(NODES{01,ikid2}([1,2]));...
         mean(NODES{01,ikid2}([3,4]))];
if (xxce(1) < xxcw(1))
  fprintf(1,'ERROR: The merge assumes "ison1" is the south box.\n')
  keyboard
end


%%% Set up the three index vectors.
J1w = [1:nedge 3*nedge+1:6*nedge];
J3w = nedge + 1:3*nedge;
J2e = 1:4*nedge;
J3e = 6*nedge:-1:4*nedge+1;

%%% Construct the solution operator.
W = (Tw(J3w,J3w) - Te(J3e,J3e))\eye(length(J3w));
U = W*[-Tw(J3w,J1w), Te(J3e,J2e)];

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

% store objects for body load contributions.
NODES{26,ibox} = W;
NODES{27,ibox} = Tw(J3w,J3w); %T_33^W
NODES{28,ibox} = Te(J3e,J3e); %T_33^e
NODES{29,ibox} = Tw(J1w,J3w);
NODES{30,ibox} = Te(J2e,J3e);
NODES{34,ibox} = [Tw(J1w,J3w); Te(J2e,J3e)]*W;
% store boundary indexing for children
NODES{31,ikid1} = J1w;
NODES{32,ikid1} = J3w;
NODES{31,ikid2} = J2e;
NODES{32,ikid2} = J3e;
% sotre indexing for exterior points
NODES{33,ibox} = indtmp;


return


function NODES = LOCAL_mergetwo_vert(NODES,ibox)

%%% Extract the relevant information from "NODES":
ikid1 = NODES{4,ibox}(1);
ikid2 = NODES{4,ibox}(2);
Ts    = NODES{23,ikid1};
Tn    = NODES{23,ikid2};
yys = NODES{13,ikid1};
yyn = NODES{13,ikid2};

nedge = length(yys)/4;

%%% Set up the three index vectors.
J1s = [1:2*nedge 3*nedge+1:4*nedge];
J3s = 2*nedge+1:3*nedge;
J2n = nedge+1:4*nedge;
J3n = nedge:-1:1;

% %%% Construct the solution operator.
W = (Ts(J3s,J3s) - Tn(J3n,J3n))\eye(length(J3n));
U = W*[-Ts(J3s,J1s), Tn(J3n,J2n)];

%%% Construct the new NfD operator;
T = [Ts(J1s,J1s), zeros(length(J1s),length(J2n));...
    zeros(length(J2n),length(J1s)), Tn(J2n,J2n)] + ...
    [Ts(J1s,J3s); Tn(J2n,J3n)] * U;

%%% Assemble the nodes in the external ring, and order them appropriately.
yyext = [yys(:,J1s),yyn(:,J2n)];
indtmp = [1:2*nedge 3*nedge+1:6*nedge 2*nedge+1:3*nedge];

% % % %%% Store away various objects 
NODES{13,ibox} = yyext(:,indtmp);    % Index vector for the external nodes.
NODES{14,ibox} = yys(:,J3s);         % Index vector for the internal nodes.
% Homogenous operators
NODES{23,ibox} = T(indtmp,indtmp);  % DtN operator   boundary to boundary
NODES{24,ibox} = U(:,indtmp);       % Solve operator boundary to int


% store objects for body load contributions.
NODES{26,ibox} = W;
NODES{27,ibox} = Ts(J3s,J3s); %T_33^s
NODES{28,ibox} = Tn(J3n,J3n); %T_33^n
NODES{29,ibox} = Ts(J1s,J3s);
NODES{30,ibox} = Tn(J2n,J3n);
NODES{34,ibox} = [Ts(J1s,J3s); Tn(J2n,J3n)]*W;



% store boundary indexing for children
NODES{31,ikid1} = J1s;
NODES{32,ikid1} = J3s;
NODES{31,ikid2} = J2n;
NODES{32,ikid2} = J3n;
% sotre indexing for exterior points
NODES{33,ibox} = indtmp;


return
