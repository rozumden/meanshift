% Copyright (c) 2017 James Pritts, Denys Rozumnyi
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
function [mu0,v,w] = test(bandwidth_type)
if nargin < 1
    bandwidth_type = 'sample';
end

rng(15762);

mu0 = [-4 4; 3  4];

c1 = 1:500;
c2 = 501:1000;

gtG = [ones(1,numel(c1)) 2*ones(1,numel(c2))];
mu = [mu0(:,gtG(c1)) mu0(:,gtG(c2))];
u = normrnd(mu,2);

switch bandwidth_type
  case 'fixed'
    meanshift = MMS.Meanshift('manifold','R2', ...
                              'bandwidth_type','fixed', ...
                              'bandwidth',3); 
  case 'sample' 
    meanshift = MMS.Meanshift('manifold','R2', ...
                              'bandwidth_type','sample', ...
                              'pilot_knn',ceil(0.3*size(u,2)));
    
end

[clust,G,modes,likelihoods,bandwidths,tracks] = meanshift.fit_and_predict(u);
meanshift.delete();

cmap0 = mat2cell(distinguishable_colors(size(mu,2))',3,ones(1,numel(gtG)));
colors0 = cmap0(gtG);

figure;
subplot(2,2,1);
hold on;
splitapply(@draw_pts,u,colors0,gtG);
hold off;
axis equal;

[G,uG] = findgroups(G);
cmap = mat2cell(distinguishable_colors(numel(uG)+1)', ...
                3,ones(1,numel(uG)+1));
colors = cmap(G);

eind1 = find(G(c1) ~= G(c1(1)));
eind2 = c2(find(G(c2) ~= G(c2(1))));

subplot(2,2,2);
hold on;
splitapply(@draw_pts,u,colors,G);
plot(u(1,eind1),u(2,eind1), ...
     'o','MarkerEdgeColor',cmap{end});
plot(u(1,eind2),u(2,eind2), ...
     'o','MarkerEdgeColor',cmap{end});
hold off; 
axis equal;

subplot(2,2,3);
hold on;
splitapply(@draw_pts,u,colors,G);
draw_circles(u(:,eind1),bandwidths(eind1));
draw_circles(u(:,eind2),bandwidths(eind2));
hold off;
axis equal;

subplot(2,2,4);
hold on;
splitapply(@draw_pts,u,colors,G);
[~,ia]=min(u(1,eind1));
plot(u(1,eind1(ia)),u(2,eind1(ia)),'rx','LineWidth',3);
plot(modes(1,eind1(ia)),modes(2,eind1(ia)),'bx','LineWidth',3);
track = tracks{eind1(ia)};
plot(track(1,:),track(2,:),'k-','LineWidth',5);
hold off;
axis equal;

keyboard;


function [] = draw_pts(u,c)
c = c';
for k = 1:size(u,2)
    plot(u(1,k),u(2,k),'Color',[c{k}]','Marker','.');
end

function [] = draw_circles(x,radii)
m = size(x,2);
t = 0:0.01:2*pi;
for k = 1:m
    A = [radii(k) 0 x(1,k); ...
         0 radii(k) x(2,k); ...
         0  0 1];
    y = A*[cos(t);sin(t);ones(1,length(t))];
    plot(y(1,:),y(2,:),'g');
end
