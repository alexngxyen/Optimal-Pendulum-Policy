% Simple value iteration, for points on a ring.
% 1 = -1 on the ring.
% Katie Byl, 6/1/20

xgoal = 0;
Vbig = 1000; % some very large initial value/penalty at each state
gamma = 1; % discount factor?  Try gamma = .99, or .9, etc.

xset = -1:.2:1;
uset = -2:2;
V = 0*xset + Vbig;
fi = find(xset == xgoal);
V(fi) = 0 % GOAL has zero value/penalty
dt = 0.05; % time step

Uopt = 0*xset; % optimal action, for each state

% First, pre-calculate all possible transitions:
Tmat = zeros(length(uset),length(xset));
for xi = 1:length(xset)
    for ui = 1:length(uset)
        xnext = xset(xi) + dt * (1 + uset(ui));
        % Make sure to "wrap" the state, if needed!
        if xnext>1
            xnext = xnext-2;
        elseif xnext<-1
            xnext = xnext+2;
        end
        % Find the closest mesh point
        [~,id] = sort(abs(xnext-xset));
        Tmat(ui,xi) = id(1); % assume it goes to nearest neighbor
    end
end

Tmat % Here is the so-called "transition matrix"

% Then, keep iterating to find the "best possible option" to perform:
for n=1:100
    Vnext = 0*V;
    for xi = 1:length(xset)
        [vi,id] = sort(V(Tmat(:,xi)));
        Vnext(xi) = gamma*vi(1) + dt*(xset(xi)~=xgoal);
        Uopt(xi) = id(1);
    end
    if abs(sum(V-Vnext)) < 10*eps
        fprintf('n = %d. No change in value!\n',n)
        break
    end
    V = Vnext
end

fprintf('The value, V, gives time-to-go, until x = xgoal. It has converged.\n')    