function x_vals = generate_xvals_v4 ( y_vals, x_offset, min_x_distance, min_y_distance, step_size )

    if (nargin < 5)
        step_size = 0.1;
    end

    %Force y_vals to be a row vector
    y_vals = y_vals(:)';
    x_vals = nan(1, length(y_vals));
    
    for i = 1:length(y_vals)
        x = x_offset;
        y = y_vals(i);
        
        y_dist = abs(y_vals - y);
        y_dist(i) = nan;
        if (any(y_dist < min_y_distance))
            idx_of_conflicting_points = find(y_dist < min_y_distance);
            associated_x_values = x_vals(idx_of_conflicting_points);
            associated_x_values(isnan(associated_x_values)) = [];
            
            if (~isempty(associated_x_values))
                xt1 = x_offset;
                xt2 = x_offset;
                done = 0;
                
                while (~done)
                    xd1 = min(abs(associated_x_values - xt1));
                    xd2 = min(abs(associated_x_values - xt2));
                    if (xd1 < min_x_distance && xd2 < min_x_distance)
                        xt1 = xt1 - step_size;
                        xt2 = xt2 + step_size;
                    else
                        done = 1;
                        if (xd1 >= min_x_distance)
                            x = xt1;
                        else
                            x = xt2;
                        end
                    end
                end
                
            end
        end
        
        x_vals(i) = x;
    end
end






































