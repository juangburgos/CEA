function mat2sql(matfile,dbfile)

try
    
% Load *.mat file in workspace
load(matfile);
% List all variables in the workspace
allmats = whos;
% Define database file name
database = dbfile;
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

%% Fill in matrices

% Open or create database file
mksqlite('open', database);
mksqlite('PRAGMA synchronous = OFF');

% Create info table
mksqlite('create table info (matname text, numrows int, numcols int)');
mksqlite('begin');
for k = 1:1:max(size(allmats))
    mksqlite(['insert into info (matname, numrows, numcols) VALUES (''', ...
        allmats(k).name, ''', ', num2str(allmats(k).size(1)), ', ', ...
        num2str(allmats(k).size(2)), ')']);
end
mksqlite('commit');

% Fill in data
for k = 1:1:max(size(allmats))
    
    table = allmats(k).name;
    eval(['mymat = ',table,';']);

    nrows = allmats(k).size(1);
    ncols = allmats(k).size(2);

    cols  = '';
    for i = 1:1:ncols
        if i < ncols
            cols = strcat(cols,['col',int2str(i),' double, ']);
        else
            cols = strcat(cols,['col',int2str(i),' double']);
        end
    end

    mksqlite(['create table ', table, ' (', cols, ')']);

    mksqlite('begin');

    cols  = '';
    for i = 1:1:ncols
        if i < ncols
            cols = strcat(cols,['col',int2str(i),', ']);
        else
            cols = strcat(cols,['col',int2str(i)]);
        end
    end
    for j = 1:1:nrows
        vals  = '';
        for i = 1:1:ncols
            if i < ncols
                vals = strcat(vals,[sprintf('%25.25f',mymat(j,i)),', ']);
            else
                vals = strcat(vals,[sprintf('%25.25f',mymat(j,i))]);
            end
        end
        mksqlite(['insert into ', table, ' (', cols, ') VALUES (',vals, ')']);
    end

    mksqlite('commit');
    
end

mksqlite('close');

% %% Check written data
% 
% mksqlite('open', database);
% 
% disp ('Lese alle Records in ein Array ein');
% 
% res = mksqlite(['SELECT * FROM ' table]);
% % sqlite> SELECT * FROM info WHERE matname = 'gamma';
% % sqlite> SELECT numrows,numcols FROM info WHERE matname = 'gamma';
% 
% mksqlite('close');

catch err
    keyboard;
end