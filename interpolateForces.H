vector ddpoint(0.0,0.0,0.0);
scalar hh = mesh.C()[1][0] - mesh.C()[0][0];

ddpoint[0] = mesh.C()[icell][0] - rr[0];
ddpoint[1] = mesh.C()[icell][1] - rr[1];
ddpoint[2] = mesh.C()[icell][2] - rr[2];

F[icell][0] = F[icell][0] + pf[0]*diracdelta(ddpoint,hh);
F[icell][1] = F[icell][1] + pf[1]*diracdelta(ddpoint,hh);
F[icell][2] = F[icell][2] + pf[2]*diracdelta(ddpoint,hh);

