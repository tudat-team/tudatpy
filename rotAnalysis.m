clc
close all
clear all

datafolder = '/home/dominic/Software/numericalAstrodynamicsTudatBundle/tudatBundle/tudat/Tudat/'
load(strcat(datafolder,'rotationFormalEstimationError.dat'));
load(strcat(datafolder,'rotationInformationMatrix.dat'));
load(strcat(datafolder,'rotationCorrelations.dat'));
load(strcat(datafolder,'rotationFormalEstimationError.dat'));
load(strcat(datafolder,'rotationInverseNormalizedCovariance.dat'));
load(strcat(datafolder,'rotationParameterNormalization.dat'));
load(strcat(datafolder,'rotationResiduals.dat'));


figure
semilogy(rotationFormalEstimationError)

figure
plot(rotationResiduals)

figure
imagesc(rotationCorrelations)

