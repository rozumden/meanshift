classdef Meanshift < handle

    properties (Access = public)
        id
        bandwidth = 'sample' % fixed/balloon/sample
        bandwidth_parameter = 0.25
        profile = 'gaussian' % gaussian/epanechnikov
        manifold = 's2'
        scalar = 'double' % double/single
        min_support = 10
        min_ratio = 1
        max_bw = 0.4
    end

    methods
        function this = Meanshift(varargin)
            if exist('cmp_argparse.m','file')
                if numel(varargin) == 1
                    varargin = KEY.class_to_struct(varargin{1});
                    varargin = [fieldnames(varargin) struct2cell(varargin)]';
                end  
                [this,~] = cmp_argparse(this,varargin{:});
            else
                warning('Cant process varargin. Helpers are not found.');
            end
            this.id = meanshift('new', 'profile',this.profile, ...
                 'manifold',this.manifold,'bandwidth',this.bandwidth, ...
                 'bandwidth_parameter', this.bandwidth_parameter, ...
                 'scalar', this.scalar, 'min_support', this.min_support, ...
                 'min_ratio', this.min_ratio,'max_bw',this.max_bw);
        end

        function [clust,lbl,modes,likelihoods] = train(this, u)
            clust = [];
            lbl = [];
            modes = [];
            likelihoods = [];
            if isempty(u)
                return;
            end
            if strcmp(this.scalar,'single')
                u = single(u);
            end
            [clust,lbl,modes,likelihoods] = meanshift('train',this.id,u);
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
            dist = meanshift('est_dist',this.id,u,clust,lbl);
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
                    text(u(1,cl),u(2,cl),u(3,cl),num2str(i),'FontSize',10,'FontWeight','bold');
                end
                outliers = lbl == 0;
                scatter3(ax1,u(1,outliers),u(2,outliers),u(3,outliers),'b.');
            end 
            plot3(modes(1,:),modes(2,:),modes(3,:),'r.','MarkerSize',15);
            for k = 1:size(modes,2)
                text(modes(1,k),modes(2,k),modes(3,k),num2str(k),'FontSize',10,'FontWeight','bold');
            end
            axis equal;
            hold off;
        end

        function delete(this)
            meanshift('delete', this.id);
        end

    end

end