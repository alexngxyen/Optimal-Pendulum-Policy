% Simple value iteration, for points on a ring.
% 1 = -1 on the ring.
% Katie Byl, 6/1/20
%
% Now, with a dt that does NOT divide nicely for the dynamics...

xgoal = 0;
Vbig = 1000; % some very large initial value/penalty at each state
gamma = 1; % discount factor?  Try gamma = .99, or .9, etc.

xset = -1:.2:1;
uset = -2:2;
V = 0*xset + Vbig;
fi = find(xset == xgoal);
V(fi) = 0 % GOAL has zero value/penalty
dt = 0.25; % time step. Try other values. 0.05, 0.01, etc...
% Make N larger as dt gets smaller, or the value will NOT converge!

Uopt = 0*xset; % optimal action, for each state

% This time, use TWO Tmats, and also save the WEIGHTS:
Tmat1 = zeros(length(uset),length(xset));
Tmat2 = Tmat1;
Wt1 = 0*Tmat1;
Wt2 = Wt1;
for xi = 1:length(xset)
    for ui = 1:length(uset)
        xnext = xset(xi) + dt * (1 + uset(ui));
        % Make sure to "wrap" the state, if needed!
        if xnext>1
            xnext = xnext-2;
        elseif xnext<-1
            xnext = xnext+2;
        end
        % Find the mesh point at or "just below" current value:
        lower_value = xnext - mod(xnext,0.2);
        upper_value = lower_value + 0.2;
        
        [~,id1] = sort(abs(lower_value-xset));
        [~,id2] = sort(abs(upper_value-xset));
        Tmat1(ui,xi) = id1(1); % assume it goes to nearest neighbor
        Tmat2(ui,xi) = id2(1); % assume it goes to nearest neighbor
        Wt2(ui,xi) = (xnext - lower_value) / 0.2;
        Wt1(ui,xi) = 1 - Wt2(ui,xi);
    end
end

Tmat1 % Now, we have N+1 = 2 Tmats, for barycentric interpolation
Tmat2

% Then, keep iterating to find the "best possible option" to perform:
for n=1:20  % Make this much LARGER, as dt gets SMALLER
    Vnext = 0*V;
    Vt = V'; % V transpose needed here, for .* operations
    for xi = 1:length(xset)
        [vi,id] = sort( Wt1(:,xi) .* Vt(Tmat1(:,xi)) + ...
            Wt2(:,xi) .* Vt(Tmat2(:,xi)) );
        
        Vnext(xi) = gamma*vi(1) + dt*(xset(xi)~=xgoal);
        Uopt(xi) = id(1);
    end
    %if abs(sum(V-Vnext)) < 10*eps,
    fprintf('n = %d. Difference in all values is %.12f\n', n,sum(abs(V-Vnext)))
    %end
    V = Vnext
end
