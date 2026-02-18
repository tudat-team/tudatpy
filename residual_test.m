clc
close all
clear all

load('/home/dominic/Downloads/prefit_modified.dat')
load('/home/dominic/Downloads/ppostfit_modified.dat')

vectorsize = length(prefit_modified);

for i=1:3
    figure(1)
    scatter(1:(vectorsize/3),log10(abs(prefit_modified(i:3:vectorsize))))
    hold on

    figure(2)
    scatter(1:(vectorsize/3),log10(abs(ppostfit_modified(i:3:vectorsize))))
    hold on
end

load('/home/dominic/Downloads/prefit_unmodified.dat')
load('/home/dominic/Downloads/ppostfit_unmodified.dat')

vectorsize = length(prefit_unmodified);

for i=1:3
    figure(3)
    scatter(1:(vectorsize/3),log10(abs(prefit_unmodified(i:3:vectorsize))))
    hold on

    figure(4)
    scatter(1:(vectorsize/3),log10(abs(ppostfit_unmodified(i:3:vectorsize))))
    hold on
end