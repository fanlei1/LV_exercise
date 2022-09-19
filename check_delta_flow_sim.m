function flow = check_delta_flow_sim(network_matrix,Part,Pven,dP_dt,G_vessels,C_vessels,P)

node_vector = unique(network_matrix(:,5:6));                               % node number (proximal and distal)
nodes = max(node_vector);
nm = zeros(nodes, 4); % 4 indicates there are maximum of 4 vessels connected to a node
for i = 1:nodes
    [row,col] = find(network_matrix(:,5:6)==i);
    nm(i,1:size(row)) = row';
end
flow.Pnode = zeros(1,nodes);
% vessels(1)=find(network_matrix(:,7)==0);
% vessels(2:length(find(network_matrix(:,9)==0))+1)=find(network_matrix(:,9)==0);
Pin = Part;
Pout = Pven;
flow.Pnode(:,1) = Pin';
for j = 2:nodes
    row = nonzeros(nm(j,:));
    if(size(row,1)==1)
        flow.Pnode(:,j)=Pout';
    elseif(size(row,1)==2)
        flow.Pnode(:,j) = (P(:,nm(j,1)).*G_vessels(:,nm(j,1))+P(:,nm(j,2)).*G_vessels(:,nm(j,2)))./ ...
            (G_vessels(:,nm(j,1))+G_vessels(:,nm(j,2)));
    elseif(size(row,1)==3)
        flow.Pnode(:,j)=(P(:,nm(j,1)).*G_vessels(:,nm(j,1))+P(:,nm(j,2)).*G_vessels(:,nm(j,2))+ ...
            P(:,nm(j,3)).*G_vessels(:,nm(j,3)))./(G_vessels(:,nm(j,1))+G_vessels(:,nm(j,2))+G_vessels(:,nm(j,3)));
    else
        flow.Pnode(:,j)=(P(:,nm(j,1)).*G_vessels(:,nm(j,1))+P(:,nm(j,2)).*G_vessels(:,nm(j,2))+ ...
            P(:,nm(j,3)).*G_vessels(:,nm(j,3))+P(:,nm(j,4)).*G_vessels(:,nm(j,4)))./ ...
            (G_vessels(:,nm(j,1))+G_vessels(:,nm(j,2))+G_vessels(:,nm(j,3))+G_vessels(:,nm(j,4)));
    end
end

flow.flow_in(:,:)= 2*(flow.Pnode(:,network_matrix(:,5))-P(:,:)).*G_vessels(:,:)*1.2;
flow.flow_out(:,:) = 2*(P(:,:)-flow.Pnode(:,network_matrix(:,6))).*G_vessels(:,:);%+C_vessels(:,:).*dP_dt(:,:);