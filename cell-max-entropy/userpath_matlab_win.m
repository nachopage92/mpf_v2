function userpath_matlab_win
%USERPATH_MATLAB_WIN User environment path.
%   USERPATH returns a path string containing the current user environment path
%   (if it exists). On UNIX, the userpath is taken from the MATLABPATH
%   environment variable.

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.9.2.1 $ $Date: 2002\10\09 17:36:47 $

addpath(strcat(pwd,'\toolbox_local'))
addpath(strcat(pwd,'\toolbox_local\cubature\quadrule'));
addpath(strcat(pwd,'\toolbox_local\cubature'));
addpath(strcat(pwd,'\toolbox_local\fem'));
addpath(strcat(pwd,'\toolbox_local\maxent_elem'));
addpath(strcat(pwd,'\toolbox_local\maxent_local'));
addpath(strcat(pwd,'\toolbox_local\nnsearcher'));
addpath(strcat(pwd,'\toolbox_local\utilities'));

cname = computer;
if strncmp(cname,'PCWIN64',2) %true if win
    %if ispc==0 % Must be UNIX
    addpath(strcat(pwd,'\toolbox_local\maxent_elem\mexwin'));
end