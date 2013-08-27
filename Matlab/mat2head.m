function mat2head(matfile,hfile)

try
    
% Load *.mat file in workspace
load(matfile);
% List all variables in the workspace
allmats = whos;
% Remove non-numerical variables
k = 1;
while 1
    if k > max(size(allmats))
        break;
    end
    if ~strcmp(allmats(k).class,'double')
        allmats(k) = [];
    else
        k = k+1;
    end
end

fileID   = fopen(hfile, 'w');

%% Fill in matrices

% Print data
for k = 1:1:max(size(allmats))
    
    
    mname = allmats(k).name;
    nrows = allmats(k).size(1);
    ncols = allmats(k).size(2);
    
    % Temporary Matlab matrix
    eval(['mymat = ',mname,';']);
    
    % Write contents to file
    fprintf(fileID, ...
     [mname,'.m = ',num2str(nrows),';\n', ...
      mname,'.n = ',num2str(ncols),';\n', ...
      mname,'.mat = (doublereal*) calloc(', ...
      mname,'.m * ',mname,'.n, sizeof(doublereal));\n', ...
      'static doublereal t_',mname,'[] = {']);

    % Print Matrix Values 
    for j = 1:1:ncols
        for i = 1:1:nrows
            if i < nrows
                fprintf(fileID,[sprintf('%25.25f',mymat(i,j)),',']);
            end
        end
        if j < ncols
            fprintf(fileID,[sprintf('%25.25f',mymat(i,j)),',\n']);
        else
            fprintf(fileID,[sprintf('%25.25f',mymat(i,j)),'};\n']);
        end
    end
    
    fprintf(fileID, ...
     [mname,'.mat = &t_',mname,'[0];\n\n\n']);

    
end

fclose(fileID);

catch err
    keyboard;
end