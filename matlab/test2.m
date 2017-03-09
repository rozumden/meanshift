addpath('~/src/meanshift/.build/Linux/x86_64/R2013a/release/cpp');

% cfg = CFG.get();
% [sqldb,imagedb] = get_dbs(get_dbs_cfg());
% cache = CASS.CidCache('12091996', imagedb);
% cache.add_dependency('meanshift',[]);
% data = cache.get('data','meanshift');
load('data');
v = data.v;

meanshift = Meanshift();
[clust,lbl,modes] = meanshift.train(v);
% p = meanshift.est_density(modes,v)
meanshift.draw(v,modes);
meanshift.draw(v,clust(:,:),lbl);
% meanshift.draw(data.v,data.modes);
