   // Let's calculate the velocity manually
   const vectorField& C = mesh.C();
   // Get current time
   const scalar t = runTime.value();
   
   // Parameters for rotation
   const scalar omega0 = 1.0;  // Base angular velocity [rad/s]
   const scalar freq = 0.1;    // Frequency of oscillation [Hz]
   scalar x0 = 0.05;     // Center of rotation x-coordinate
   scalar y0 = 0.05;     // Center of rotation y-coordinate

    // Time-dependent angular velocity
    // Example: Sinusoidally varying rotation speed
    const scalar omega = omega0 * (1.0 + Foam::sin(2.0*3.14*freq*t));

   forAll(U, i)
   {
       const scalar x = C[i].x();
       const scalar y = C[i].y();
       
       // Calculate rotational velocity components
       U[i].x() = -omega * (y - y0);
       U[i].y() = omega * (x - x0);
       U[i].z() = 0.0;
   }
