%% select compiler
% mex -setup C++

%% compile mex file
clc
mex -g HiLoNLM_mex_new_v2.cpp  %updated for seperate denoise and despeckle

