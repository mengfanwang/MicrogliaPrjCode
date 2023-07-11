% clc;clear;close all;
% load status2.mat

function plotHull(line_ind,face,vertex)
figure(1);
close(1);
figure(1);
plot3(line_ind(1,:),line_ind(2,:),line_ind(3,:),'.k','MarkerSize',10);
hold on;
plot3(0,0,0,'.r','MarkerSize',10);
% plot3(dist_minvertex(1,:),dist_minvertex(2,:),dist_minvertex(3,:),'.r','MarkerSize',10);
% plot3(vertex(1,:),vertex(2,:),vertex(3,:),'.g','MarkerSize',10);
for ff = 1:length(face)
    if ~isempty(face{ff})
        trisurf([1 2 3],face{ff}.vertex(1,:),face{ff}.vertex(2,:),face{ff}.vertex(3,:));
    end
end
