# SN_RKHS

This file contains the Matlab code for the paper "AN RKHS APPROACH FOR PIVOTAL INFERENCE IN FUNCTIONAL LINEAR REGRESSION" by Holger Dette and Jiajun Tang

This file contians:

(1) sample_code.m: A sample simulation code for computing empirical rejection probabilities. This simulation takes several hours to finish.

(2) TandV.m: a function that computes \hat{T}_n and \hat{V}_n in the paper.

(3) bike.m: code for real data example Bike-sharing.

(4) bikeshare.txt: The Bike-sharing data of Captial Bike Sharing (CBS), obtained from R package ISLR2. 

References of this data set:

Fanaee-T, H. and Gama, J. (2014). Event labeling combining ensemble detectors and background knowledge. Progress in Artificial Intelligence, 2, 113–127.

James, G., Witten, D., Hastie, T. and Tibshirani, R. (2021). An Introduction to Statistical Learning: With Applications in R. New York: springer.

The code uses Matlab package Chebfun available at https://www.chebfun.org/. Please download this package and put it in the same folder as the above files and set this folder as the working directory.
