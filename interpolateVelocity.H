vector ddpoint(0.0,0.0,0.0);
scalar hh = mesh.C()[1][0] - mesh.C()[0][0];

ddpoint[0] = mesh.C()[icell][0] - rr[0];
ddpoint[1] = mesh.C()[icell][1] - rr[1];
ddpoint[2] = mesh.C()[icell][2] - rr[2];

pu[0] = pu[0] + U[icell][0]*diracdelta(ddpoint,hh)*(hh*hh*hh);
pu[1] = pu[1] + U[icell][1]*diracdelta(ddpoint,hh)*(hh*hh*hh);
pu[2] = pu[2] + U[icell][2]*diracdelta(ddpoint,hh)*(hh*hh*hh);

