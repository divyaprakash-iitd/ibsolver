// Create the cloud
Cloud<passiveParticle> ibpoints(mesh, "pointsCloud", false);

// Add particles
vector pos1(0.025, 0.075, 0.005);
ibpoints.addParticle(new passiveParticle(mesh, pos1));

