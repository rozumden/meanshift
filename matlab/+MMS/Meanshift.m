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
classdef Meanshift < handle
    properties (Access = public)
        % floating point precision 
        scalar = 'double' % double/single

        manifold = 'R2'
        
        profile = 'gaussian' % gaussian/epanechnikov
 
        % mode clustering parameters
        min_support = 21
        min_ratio = 1

        % bandwidth estimation parameters 
        bandwidth_type = 'fixed' % fixed/balloon/sample
        bandwidth = 0.2
        max_bandwidth = 0.3
        pilot_knn = 5
    end

    properties (Access = private)
        id
    end

    methods
        function this = Meanshift(varargin)
            addpath('~/opt/mex');
            if exist('+MMS/cmp_argparse.m','file')
                if numel(varargin) == 1
                    varargin = KEY.class_to_struct(varargin{1});
                    varargin = [fieldnames(varargin) struct2cell(varargin)]';
                end  
                [this,~] = MMS.cmp_argparse(this,varargin{:});
            else
                warning('Cant process varargin. Helpers are not found.');
            end
          
            this.id = meanshift('new', ...
                                'scalar', this.scalar, ...
                                'manifold', this.manifold, ...
                                'profile', this.profile, ...
                                'min_support', uint32(this.min_support), ...
                                'bandwidth_type', this.bandwidth_type, ...
                                'bandwidth', this.bandwidth, ...
                                'max_bandwidth', ...
                                this.max_bandwidth, ...
                                'pilot_knn',uint32(this.pilot_knn));
        end

        function [clust,labeling,modes,likelihoods,bandwidths,tracks] = ...
                fit_and_predict(this,u)
            if isempty(u)
                clust = [];
                labeling = [];
                modes = [];
                likelihoods = [];
                tracks = [];

                return;
            end
            u = feval(this.scalar,u);
            [clust,labeling,modes,likelihoods,bandwidths,tracks] = ...
                meanshift('fit_and_predict',this.id,u);
            labeling = double(labeling);
            labeling(labeling == 0) = nan;
            labeling = findgroups(labeling);
        end

        function p = est_density(this,modes,u)
            if isempty(modes) | isempty(u)
                return;
            end
            p = meanshift('est_density',this.id,modes,u);
        end

        function dist = est_dist(this,u,clust,lbl)
            dist = [];
            if isempty(clust) | isempty(u)
                return;
            end
            
            dist = meanshift('est_dist',this.id,u, ...
                             uint32(clust),uint32(lbl));
        end

        function draw(this,u,modes,lbl,ax1)
            if strcmp(this.scalar,'single')
                modes = double(modes);
            end
            if nargin < 4 
                lbl = [];
                figure;
                ax1 = gca;
            end
            if nargin == 4
                figure;
                ax1 = gca;
            end
            if isempty(modes)
                warning('No modes found!');
                return;
            end
            hold on;
            [xx,yy,zz] = sphere(120);
            w = [xx(:) yy(:) zz(:)]';
            p = this.est_density(w,u);
            p2 = reshape(p,size(xx));
            C = ones(size(p2,1),size(p2,2),3);
            C(:,:,1) = 1-((p2)-min(p))/(max(p)-min(p));
            % C = ones(size(xx,1),size(xx,2),3);
            surf(ax1,xx,yy,zz,C,'EdgeColor','None');
            axis equal; 
            cameratoolbar;
            if isempty(lbl)
                scatter3(ax1,u(1,:),u(2,:),u(3,:),'b.');
            else
                ulbl = setdiff(unique(lbl),0);
                mpdc = distinguishable_colors(numel(ulbl)+1);
                for i = 1:numel(ulbl)
                    cl = lbl == ulbl(i);
                    scatter3(ax1,u(1,cl),u(2,cl),u(3,cl),'CData',mpdc(i+1,:));
                    text(u(1,cl),u(2,cl),u(3,cl),num2str(i),'FontSize',100,'FontWeight','bold');
                end
                outliers = lbl == 0;
                scatter3(ax1,u(1,outliers),u(2,outliers),u(3,outliers),'b.');
            end 
            plot3(modes(1,:),modes(2,:),modes(3,:),'r.','MarkerSize',15);
            for k = 1:size(modes,2)
                text(modes(1,k),modes(2,k),modes(3,k),num2str(k),'FontSize',100,'FontWeight','bold');
            end
            axis equal;
            hold off;
        end

        function delete(this)
            meanshift('delete', this.id);
        end

    end

end