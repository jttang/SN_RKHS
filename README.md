# SN_RKHS

This package contains the Matlab code for the paper "AN RKHS APPROACH FOR PIVOTAL INFERENCE IN FUNCTIONAL LINEAR REGRESSION" by Holger Dette and Jiajun Tang

The files are:

(1) sample_code.m: A sample simulation code for computing empirical rejection probabilities.

(2) TandV.m: a function that computes \hat{T}_n and \hat{V}_n in the paper.

(3) bike.m: code for real data example Bike-sharing.

(4) bikeshare.txt: The Bike-sharing data of Captial Bike Sharing (CBS), obtained from R package ISLR2. 

References for the Bike-sharing data set:

Fanaee-T, H. and Gama, J. (2014). Event labeling combining ensemble detectors and background knowledge. Progress in Artificial Intelligence, 2, 113–127.

James, G., Witten, D., Hastie, T. and Tibshirani, R. (2021). An Introduction to Statistical Learning: With Applications in R. New York: springer.

The code uses Matlab package Chebfun available at https://www.chebfun.org/. Please download and source this package before running using the code.
