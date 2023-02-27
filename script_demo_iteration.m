%% script_demo_iteration
%
% This script runs one iteration of the "simplified MOCBA" algorithm,
% i.e., computes the proportions for a given set of means and variances.
%
% Running this script should produce the following output:
%  |
%  | >> script_demo_iteration
%  |
%  |  alpha =
%  |
%  |      0.2667    0.2000    0.4000    0.1333
%  |

input_set_size = 4;
nb_objectives = 2;

obj_mean = [0.0 0.1;
            0.1 1.0;
            1.0 0.2;
            1.1 0.3];

obj_var = [1 2;
           3 4;
           5 6;
           7 8] / 20;

assert (isequal (size (obj_mean), [input_set_size nb_objectives]));
assert (isequal (size (obj_var),  [input_set_size nb_objectives]));

alpha = sMOCBA_iteration (obj_mean, obj_var)  %#ok<NOPTS> 
