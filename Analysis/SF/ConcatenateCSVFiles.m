%clc;clear;close all;

%[filename,pathname] = uigetfile('*.*', 'Load file', 'MultiSelect', 'on')

A1 = importdata('StepSameInd-SKS007.csv')
A2 = importdata('StepSameInd-SKS008.csv')

%A= importdata(filename(1,1))