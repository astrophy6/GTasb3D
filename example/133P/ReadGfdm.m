%% read binary file from GTasb2D simulation
function [time, f, data, data_gas] = ReadGfdm(filename, num_node, flag_gas, node_data)

fileID = fopen(filename,'r');

time = fread(fileID,1,'double'); % time
f = fread(fileID,1,'double'); % true anomaly
data = fread(fileID,[3,num_node],'double'); 
data = data';

fclose(fileID);

if flag_gas
    fileID = fopen(filename+".gas",'r');
    data_gas = zeros([num_node,4]);
    for i = 1:num_node
        data_gas(i,1) = fread(fileID,1,'double');
        if node_data(i,1) < 2
            data_gas(i,2) = fread(fileID,1,'double');
        else
            data_gas(i,2:4) = fread(fileID,[1,3],'double');
        end
    end
    fclose(fileID);
else
    data_gas = 0;
end

end



