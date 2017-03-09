function [mu0,v,w] = test(bw,min_card,min_dist)

if nargin == 0
	bw = 0.2;
    min_card = 10;
    min_dist = 0.1;
end

addpath('~/src/meanshift/.build/Linux/x86_64/R2013a/release/cpp');

mu0 = [-10  0; ...
       3  5];
c1 = 1:1000;
c2 = 1001:2000;
mu = [mu0(:,ones(1,numel(c1))) mu0(:,2*ones(1,numel(c2)))];
u = normrnd(mu,1);
[v,w] = MEANSHIFT.meanshift(u,bw,min_card,min_dist);

figure;
subplot(1,3,1);
hold on;
plot(u(1,c1),u(2,c1),'b.');
plot(u(1,c2),u(2,c2),'r.');
hold off;
axis equal;


subplot(1,3,2);
hold on;
plot(u(1,c1),u(2,c1),'b.');
plot(u(1,c2),u(2,c2),'r.');
plot(v(1,:),v(2,:),'g.', ...
     'MarkerSize',15);
hold off;
axis equal;

subplot(1,3,3);
hold on;
plot(u(1,c1),u(2,c1),'b.');
plot(u(1,c2),u(2,c2),'r.');
plot(w(1,:),w(2,:),'g.', ...
     'MarkerSize',15);
hold off;
axis equal;