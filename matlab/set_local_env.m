
addpath('./data_utils/');
addpath('./gen_utils/');
addpath('./out_utils/');

if contains(computer,'mac','IgnoreCase',true)
    % On my local machine
    data_dir = '/Users/stharakan/Documents/local_data/';
    libsvm_dir = '/Users/stharakan/Documents/lib/libsvm-master/matlab';
else
    % on a tacc computer
    data_dir = '/work/03158/tharakan/data/';
    libsvm_dir = '/work/03158/tharakan/lib/libsvm-master/matlab';
end
addpath(libsvm_dir);

runfile_dir='./../runfiles/';

