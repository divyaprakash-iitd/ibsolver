   // Let's calculate the velocity manually
   const vectorField& C = mesh.C();
   
   // Parameters for rotation
   scalar omega = 1.0;  // Angular velocity [rad/s]
   scalar x0 = 0.05;     // Center of rotation x-coordinate
   scalar y0 = 0.05;     // Center of rotation y-coordinate
   
   forAll(U, i)
   {
       const scalar x = C[i].x();
       const scalar y = C[i].y();
       
       // Calculate rotational velocity components
       U[i].x() = -omega * (y - y0);
       U[i].y() = omega * (x - x0);
       U[i].z() = 0.0;
   }