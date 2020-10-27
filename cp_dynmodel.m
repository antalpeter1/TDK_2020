function qd = cp_dynmodel(q,u,par,fr)
phi = q(1);
phid = q(2);
x = q(3);
xd = q(4);
m = par(1); m_0 = par(2); k_m = par(3); R = par(4); L_a = par(5);
R_a = par(6); k_e = par(7); g = par(8); s = par(9);

qd(1,1) = phid; 

qd(2,1) = ((m+m_0)*g*sin(phi) - (k_m/R/R_a*(u-k_e/R*xd)-fr(1)*...
    tanh(xd*fr(2)))*cos(phi) -...
     m*s*phid^2*sin(phi)*cos(phi)) / (4/3*s*(m+m_0) - m*s*cos(phi)^2);

qd(3,1) = xd; 
 
qd(4,1) = (k_m/R/R_a*(u-k_e/R*xd)-fr(1)*tanh(xd*fr(2))- 3/4*m*g*sin(phi)*cos(phi) +...
    m*s*phid^2*sin(phi)) / (m+m_0-3/4*m*cos(phi)^2);