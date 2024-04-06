function s = stats(omega, outputUnits)
arguments
    omega (:,2) double
    outputUnits (1,1) string {mustBeMember(outputUnits, ["deg", "rad"])} = "deg"
end
% Author: Tom Shlomo, ACLab BGU, 2020

G = size(omega,1);
m = zeros(G,1);
for i=1:G
    d = angle_between(omega(i,:), omega);
    d(i) = [];
    m(i) = min(d);
end
s.max_min = max(m);
s.avg_min = mean(m);

if outputUnits == "deg"
    s.max_min = s.max_min*180/pi;
    s.avg_min = s.avg_min*180/pi;
end

end

