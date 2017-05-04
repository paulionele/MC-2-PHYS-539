function pairElecDir = testMC(theta,phi,eta,chi)

a = [1,theta,phi]; % r theta phi


modDir = [sin(chi)*cos(eta),sin(chi)*sin(eta),cos(chi)]';

transDir1 = [sin(a(2))*cos(a(3)),   cos(a(2))*cos(a(3)),    -sin(a(3))];
transDir2 = [sin(a(2))*sin(a(3)),   cos(a(2))*sin(a(3)),     cos(a(3))];
transDir3 = [cos(a(2)),             -sin(a(2)),                      0];

transDir = [transDir1; transDir2; transDir3];
newDir = transDir*modDir;

sphereRho = (newDir(1)^2 + newDir(2)^2 + newDir(3)^2)^.5;
sphereTheta = acos(newDir(3)/sphereRho);
spherePhi = atan2(newDir(2),newDir(1)); 
%atan2 give four-quadrant arctangent, from zero to pi then -pi to zero
%meaning pi and -pi are located at the same position on the unit circle.
%We want values from 0 to 2pi for phi, but we need to maintain the proper
%postioning of the negative returned values.  We add 2pi to any negative
%values and that will complete the circle from 0 to 2pi

if abs(sphereTheta) < 1E-5 %numerical solutions of pi are not exactly zero, round down.
    sphereTheta = 0;
end
if abs(spherePhi) < 1E-5
    spherePhi = 0;
end

%We want values from 0 to 2pi for phi, but we need to maintain the proper
%postioning of the negative returned values.  We add 2pi to any negative
%values and that will complete the circle from 0 to 2pi
if spherePhi < 0
    spherePhi = spherePhi + 2*pi; 
end

pairElecDir = [sphereTheta; spherePhi;];