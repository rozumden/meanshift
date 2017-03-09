function [] = meanshift_init(opt_path)
if nargin < 1
    opt_path = '~/opt/';
end

[cur_path, name, ext] = fileparts(mfilename('fullpath'));

addpath(cur_path);

if ~exist('meanshift.mexa64','file')
    addpath([opt_path 'meanshift/matlab']);
end