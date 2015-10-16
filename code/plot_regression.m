function [ success ] = plot_regression(in,out)
%PLOT_RESULTS Plotting of the GPExp simulation results
%   Required inputs:
%   in : the input structure of GPExp
%   out: the output (results) structure of GPExp

%% 3D plot of the results:
x=in.x;
y=in.y;
z = out.ndgrid.z;
y_gp = out.ndgrid.y_gp;
s2 = out.ndgrid.s2;


nv = size(in.considered_inputs,1);

z_shaped = zeros(in.Ngrid,nv);
for i = 1:nv
    if nv > 1
        yp_shaped=reshape(y_gp,ones(1,nv)*in.Ngrid);   % make a n-D matrix
        z_shaped(:,i)=linspace(min(z(:,i)),max(z(:,i)),in.Ngrid);
    else
        z_shaped = z;
    end
end

if nv==1
    f = [y_gp+2*sqrt(s2); flipdim(y_gp-2*sqrt(s2),1)];
    fill([z; flipdim(z,1)], f, [7 7 7]/8)
    hold on; plot(z, y_gp); plot(x, y, '+')
    xlabel(in.considered_inputs); ylabel(in.considered_output);
elseif nv==2
    hold on;
    surf(z_shaped(:,1),z_shaped(:,2),yp_shaped);
    plot3(x(:,1),x(:,2), y, '+')
    xlabel(in.considered_inputs(1)); ylabel(in.considered_inputs(2)) ; zlabel(in.considered_output);
    grid on
    view(45,25);
elseif nv > 2
    if isfield(in,'idx_xy') %if the x and y axes are imposed
        idx_sort = [in.idx_xy, 1:(min(in.idx_xy)-1),(min(in.idx_xy)+1):(max(in.idx_xy)-1),(max(in.idx_xy)+1:nv)] ;
    else % Sort the variables in terms of their weights or lengthscale in absolute value:
        if ~isempty(strfind(in.covfunction{:},'ard'))
            [aa idx_sort] = sort(out.hypcov(1:end-1));
        else
            [aa idx_sort] = sort(abs(out.weights),'descend');
        end
    end
    % Set the nv-2 least relevant variables yp_shaped their median value:
    med = median(z);
    % Take a slice of the hypercube yp_shaped plot the two main relevant
    % variables:
    y_surf = permute(yp_shaped,idx_sort);
    med=med(idx_sort);
    itp=[];
    plottext = '';
    for i=1:nv
        if i==1 || i==2
            itp=[itp,',:'];
        else
            indi=idx_sort(i);
            vec = z_shaped(:,indi);
            [dd,indt] = min(abs(vec-med(i)));
            itp=[itp,[',',num2str(indt)]];
            plottext = strcat(plottext,{' '},in.considered_inputs(indi), {' = '}, num2str(vec(indt)), {'; '});
        end
    end
    eval([ 'y_surf = y_surf(' itp(2:end) ');'])
    vecx = z_shaped(:,idx_sort(1));
    vecy = z_shaped(:,idx_sort(2));
    surf(vecy,vecx,y_surf);
    ylabel(in.considered_inputs(idx_sort(1))); xlabel(in.considered_inputs(idx_sort(2))) ; zlabel(in.considered_output); title(plottext);
    grid on
    view(45,25);
end
success = true;
end