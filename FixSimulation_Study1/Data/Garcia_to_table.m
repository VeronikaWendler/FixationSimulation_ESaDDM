% Load the MAT file
load('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/Data/Garcia_Eye_for_Simulation.mat');

% data = TB.Garcia;
% 
% fields = fieldnames(data);
% 
% GarciaTable = table();
% 
% for i = 1:length(fields)
%     fieldName = fields{i};
%     fieldData = data.(fieldName);
%     
%     if isvector(fieldData) && size(fieldData, 1) == 1
%         fieldData = fieldData'; % Transpose to column vector
%     end
%     
%     if iscell(fieldData)
%         fieldData = string(fieldData);
%     end
%     
%     GarciaTable.(fieldName) = fieldData;
% end
% 
% TB.Garcia = GarciaTable;
% 
% save('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/Data/GarciaData_Transformed.mat', 'TB');



data = TB.Garcia;
fields = fieldnames(data);
GarciaTable = table();

for i = 1:length(fields)
    fieldName = fields{i};
    fieldData = data.(fieldName);
    if isvector(fieldData) && size(fieldData, 1) == 1
        fieldData = fieldData'; % Transpose 
    end

    if iscell(fieldData)
        fieldData = string(fieldData);
    end
    if (strcmp(fieldName, 'p1') || strcmp(fieldName, 'p2')) && isnumeric(fieldData)
        fieldData = fieldData / 100;  % Convert
    end

    GarciaTable.(fieldName) = fieldData;
end

TB.Garcia = GarciaTable;
save('D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/Analysis_Simulation_replication/simulation/Data/GarciaData_Transformed.mat', 'TB');