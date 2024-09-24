
nodes = [2 3 0;0 2 1];
disp = [0 0.15 -0.1;0 0.2 0.1];
% Assignment 7 Q2

x_dom = 0:0.01:3;
y_dom = 0:0.01:2;



z_x = zeros(length(y_dom),length(x_dom));
z_y = zeros(length(y_dom),length(x_dom));

for i = 1:length(y_dom)
    for j = 1:length(x_dom)
        z = inst_disp(x_dom(j),y_dom(i));
        z_x(i,j) = z(1);
        z_y(i,j) = z(2);
    end
end


[X,Y] = meshgrid(x_dom,y_dom);
Z = sqrt(z_x.^2 + z_y.^2);
contourf(X,Y,Z,10,'-.','ShowText','on')
xlabel('X-Axis')
ylabel('Y-Axis')
colorbar

A = 0.5*det(transpose([1 1 1;nodes]));
H = (0.5/A)*[nodes(2,2)-nodes(2,3) 0 nodes(2,3)-nodes(2,1) 0 nodes(2,1)-nodes(2,2) 0; 
    0 nodes(1,3)-nodes(1,2) 0 nodes(1,1)-nodes(1,3) 0 nodes(1,2)-nodes(1,1); 
    nodes(1,3)-nodes(1,2) nodes(2,2)-nodes(2,3) nodes(1,1)-nodes(1,3) nodes(2,3)-nodes(2,1) nodes(1,2)-nodes(1,1) 
    nodes(2,1)-nodes(2,2)];
node_disp = transpose([disp(1,1) disp(2,1) disp(1,2) disp(2,2) disp(1,3) disp(2,3)]);
elementstrain = H*node_disp


function ele_disp = inst_disp(x,y)

        nodes = [2 3 0;0 2 1];
        disp = [0 0.15 -0.1;0 0.2 0.1];
        
        A = 0.5*det(transpose([1 1 1;nodes]));

        N1 = (0.5/A)*det(transpose([1 1 1;x nodes(1,2:3);y nodes(2,2:3)]));
        N2 = (0.5/A)*det(transpose([1 1 1;nodes(1,1) x nodes(1,3);nodes(2,1) y nodes(2,3)]));
        N3 = (0.5/A)*det(transpose([1 1 1;nodes(1,1:2) x;nodes(2,1:2) y]));
        
        N = [N1 0 N2 0 N3 0;0 N1 0 N2 0 N3];
        d = transpose([disp(1,1) disp(2,1) disp(1,2) disp(2,2) disp(1,3) disp(2,3)]);
        
        if N1 >=0 && N2 >=0 && N3 >=0
            ele_disp = N*d;
        else 
            ele_disp = [NaN NaN];
        end
end


