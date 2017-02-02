# SPLS
Stochastic Optimization for Multiview Representation Learning using Partial Least Squares
To cite this code:
@inproceedings{arora2016stochastic,
  title={Stochastic optimization for multiview representation learning using partial least squares},
  author={Arora, Raman and Mianjy, Poorya and Marinov, Teodor},
  booktitle={Proceedings of The 33rd International Conference on Machine Learning},
  pages={1786--1794},
  year={2016}
}

To run the code, please make the following folders: "CODE", "DATA", "PAGE", "PLOTS", "REPORT", 
and download the code to the "CODE" folder.

The actual data should be in "DATA" folder, in the form of a struct with the following filds:
   data.view1.training, 
   data.view1.tuning, 
   data.view1.testing, 
   data.view2.training, 
   data.view2.tuning, 
   data.view2.testing
each of which is a "di x Nj" matrix, i={1,2}, j={1,2,3}, where Nj is the number of samples and dj is the dimension.
In addition to the data, a "perm" matrix should be included in the "DATA" folder which basically stores random 
permutations over the samples. For example, this is a code to generate such a permutation matrix:
   perm=zeros(1000,N);
   for i=1:1000
     perm(i,:)=randperm(N);
   end
   save permdata.mat perm

Please see the "Demo.m" for a simple demonstration.
