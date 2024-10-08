x_domain=0:0.01:3;
y_domain=0:0.01:2;

z_x= zeros(length(y_domain),length(x_domain));
z_y= zeros(length(y_domain),length(x_domain));

for i=1:length(y_domain)
    for j=1:length(x_domain)
        z= inst_disp(x_domain(j),y_domain(i));
        z_x(i,j)=z(1);
        z_y(i,j)=z(2);
    end
end

total_A=0.5*det(transpose([1 1 1;nodes]));
H= (0.5/total_A)*[nodes(2,2)-nodes(2,3) 0 nodes(2,3)-nodes(2,1) 0 nodes(2,1)-nodes(2,2) 0; 
    0 nodes(1,3)-nodes(1,2) 0 nodes(1,1)-nodes(1,3) 0 nodes(1,2)-nodes(1,1); 
    nodes(1,3) nodes(1,2) nodes(2,2)-nodes(2,3) nodes(1,1)-nodes(1,3) nodes(2,3)-nodes(2,1) nodes(1,2)-nodes(1,1) nodes(2,1)-nodes(2,2)];
node_disp=transpose([disps(1,1) disps(2,1) disps(1,2) disps(2,2) disps(1,3) disps(2,3)]);
ele_strain = H*node_disp

function ele_disp = inst_disp(x,y)
nodes = [2 3 0; 0 2 1];
disps= [0 0.1 -0.1; 0.1 0.2 0.15];

A= 0.5*det(transpose([1 1 1; nodes]));

N1 = (0.5/A)*det(transpose([1 1 1;x nodes(1,2:3);y nodes(2,2:3)]));
N2 = (0.5/A)*det(transpose([1 1 1;nodes(1,1) x nodes(1,3);nodes(2,1) y nodes(2,3)]));
N3 = (0.5/A)*det(transpose([1 1 1;nodes(1,1:2) x;nodes(2,1:2) y]));
 
 N = [N1 0 N2 0 N3 0;0 N1 0 N2 0 N3];
 d = transpose([disps(1,1) disps(2,1) disps(1,2) disps(2,2) disps(1,3) disps(2,3)]);
 
 if N1 >=0 && N2 >=0 && N3 >=0
 ele_disp = N*d;
 else 
 ele_disp = [NaN NaN];
 end
end

