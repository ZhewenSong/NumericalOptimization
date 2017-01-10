tmax = 20;
dt = 0.01;
N = tmax/dt;
traj = zeros(4,2);
traj(1,:)=[1,1];
traj(2,:)=[-1,1];
traj(3,:)=[-1,-1];
traj(4,:)=[1,-1];
Traj = zeros(N,4,2);

clockwise = [1, 1.3, -1, -1.2];
for t=1:N
    vel = zeros(4,2); % initialize velocity for every time step!!!
    for i=1:4
        for j=1:4
            if i~=j
                % calculating the difference of the two points
                delta = traj(j,:) - traj(i,:); 
                % adding velocity vector for each neighbor point
                vel(i,:) = vel(i,:) + [delta(2), -delta(1)] * clockwise(j) / norm(delta)^2;
            end
        end
    end
    for i=1:4
        % update trajectory
        Traj(t,i,:) = traj(i,:);
        traj(i,:) = traj(i,:) + vel(i,:)*dt;  
    end
end
for i=1:4
    plot(Traj(:,i,1),Traj(:,i,2),'o');
    hold on
end
axis equal;
            