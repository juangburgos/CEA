clear all
close all
clc

load('QU_Controller_Parameters.mat');
database = 'controller_matrices.db';

%% Fill in matrices

table = 'testmat';
mymat = double(Aj);

mksqlite('open', database);
mksqlite('PRAGMA synchronous = OFF');

nrows = size(mymat,1);
ncols = size(mymat,2);

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

mksqlite('close');

%% Check written data

mksqlite('open', database);

disp ('Lese alle Records in ein Array ein');

res = mksqlite(['SELECT * FROM ' table]);

% Datenbank weder schliessen
mksqlite('close');