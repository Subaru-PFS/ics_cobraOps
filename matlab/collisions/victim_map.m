function victim_map(collision_output,title_prefix,title_suffix)
% usage: victim_map(simfun_output.Coll,title_prefix=[],title_suffix=[])
%
% generates a collision victim map from result of detectCollisionsSparse
% 
% for convenience, this auto-detects the Coll structure from simFun output.

    % first see if this is simFun output.  If so, find the detectCollisionsSparse output
    if isfield(collision_output,'Coll')
        Coll = collision_output.Coll;
    else
        Coll = collision_output;
    end
    
    if ~exist('title_prefix','var'), title_prefix = []; end;
    if ~exist('title_suffix','var'), title_suffix = []; end;
    figure
    % get number of trajectory steps in simulation
    N_traj_steps = size(Coll.type,2);

    % vector of ID's {1:2394}
    colliders = find(Coll.V > 0);
    % vector of NN ID's {1:13782}
    collidersNN = nonzeros(Coll.rcindx(colliders,colliders));

    % use a four-color jet colormap
    colormap(jet(4));
    % list out the collision types corresponding to the colors defined
    cb_TickLabels={'none','elbow','arm','fiber'};

    % generate map.  need to increment by 1 to escape matlab's 0-1 degeneracy
    image(Coll.type(collidersNN,:).'+1);
    title([title_prefix 'Victimization map -- unvictimized perpretrators show only ''none''' ...
           title_suffix], 'interpreter','none');
    ylabel('\leftarrow trajectory step');
    xlabel('cobra ID');
    grid on;

    % by default, only look at the last 10 steps
    ylim([-10.5 0.5]+N_traj_steps);

    % set up the colorbar
    cb = lcolorbar(cb_TickLabels, ...
                   'Location','vertical',...
                   'TitleString','victimization type');
    cb.YTickLabelRotation = 90;
    cb.Position           = [.95 .15 .025 .77];

    % fix the x axis ticklabels and position -- must come after colorbar
    ax = gca;
    ax.XTick              = 1:length(collidersNN);
    ax.XTickLabel         = cellstr(num2str(sort(Coll.row(collidersNN))));
    ax.XTickLabelRotation = 90;
    ax.Position           = [.03 .15 .9 .77];
