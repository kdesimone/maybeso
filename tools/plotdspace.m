function plotdspace( dspace, param )

% PLOTDSPACE  Plot proxy decision space
%
%     usage:  plotdspace( dspace, param )
%
%     input arguments
%         'dspace' is a struct that describes the proxy decision space; it is the return argument from calcdspace.m
%         'param' is a 1 x 5 matrix of GDD function parameters; it is the return argument from fitgdd.m

% set default arguments
if nargin<2, param = []; end

% don't plot cells with too few trials
dspace.pmat(dspace.nmat<10) = NaN;

% trim response matrix (remove elements with no trials)
mi = min(isnan(dspace.pmat),[],1);
mj = min(isnan(dspace.pmat),[],2);
i1 = min( find(mi==0,1,'first'), find(mj==0,1,'first') );
i2 = max( find(mi==0,1,'last'),  find(mj==0,1,'last')  );
dspace.dlist = dspace.dlist(i1:i2);
dspace.pmat = dspace.pmat(i1:i2,i1:i2);
dspace.kmat = dspace.kmat(i1:i2,i1:i2);
dspace.nmat = dspace.nmat(i1:i2,i1:i2);

% make formatted plot of response matrix
imagesc( dspace.dlist, dspace.dlist, dspace.pmat, [ 0 1 ]) ; hold on;
nancolor([ 1 1 1 ]);
set(gca,'FontName','helvetica','FontWeight','bold','FontSize',24);
xlabel 'decision variable 1'
ylabel 'decision variable 2'
axis xy square

if ~isempty(param)
    
    % get axis limits
    xlim = get(gca,'XLim');
    ylim = get(gca,'YLim');
    hold on;
    
    % decode parameters
    theta1 = param(1);
    theta2 = param(2);
    delta1 = param(3);
    delta2 = param(4);
    
    % plot decision line 1:  [ x y ]*[ cos(theta1) sin(theta1) ]' = delta1; see Pritchett and Murray's (2015) equation (2)
    if abs(sin(theta1))>0.001
        h(1) = plot(xlim,(delta1-xlim*cos(theta1))/sin(theta1),'r-');
    else
        % special case:  vertical line
        h(1) = plot([ delta1 delta1 ]*sign(cos(theta1)),ylim,'r-');
    end
    
    % plot decision line 2:  [ x y ]*[ cos(theta2) sin(theta2) ]' = delta2; see Pritchett and Murray's (2015) equation (2)
    if abs(sin(theta2))>0.001
        h(2) = plot(xlim,(delta2-xlim*cos(theta2))/sin(theta2),'r-');
    else
        % special case:  vertical line
        h(2) = plot([ delta2 delta2 ]*sign(cos(theta2)),ylim,'r-');
    end
    
    % adjust appearance
    set(h,'LineWidth',3);
    axis([ xlim ylim ]);
    hold off;
    
end

end


function nancolor( rgb )

% NANCOLOR  Set color of NaN elements in current image
%
%     usage:  nancolor( rgb )
%
%     input argument
%         'rgb' is a 1 x 3 matrix that specifies the color of NaN elements

% set default color
if nargin<1, rgb = [ 1 1 1 ]; end

% get image data
h = get(gca,'Children');
h = findobj(h,'Type','image');
im = get(h,'CData');

% get colormap
map = colormap;
n = size(map,1);

% set first row of colormap
map(1,:) = rgb;
colormap(map);

% get data range and step size
mn = min(im(:));
mx = max(im(:));
step = (mx-mn)/(n-2);

% set color range
set(gca,'CLim',[ mn-step mx ]);

end
