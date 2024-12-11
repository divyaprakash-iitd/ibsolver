
scalar diracdelta(vector x, scalar h) {
    scalar dirac = 1.0;
    scalar dd;

    for (int ii = 0; ii < 3; ++ii) {
        scalar r = x[ii] / h;

        scalar phi;
        if (Foam::mag(r) <= 1.0) {
            phi = (3.0 - 2.0 * Foam::mag(r) + Foam::sqrt(1.0 + 4.0 * Foam::mag(r) - 4.0 * r * r)) / 8.0;
        } else if (mag(r) > 1.0 && mag(r) <= 2.0) {
            phi = (5.0 - 2.0 * Foam::mag(r) - Foam::sqrt(-7.0 + 12.0 * Foam::mag(r) - 4.0 * r * r)) / 8.0;
        } else {
            phi = 0.0;
        }

        // if (Foam::mag(r) <= 2) {
        //     phi = 0.25 * (1+Foam::cos(Foam::constant::mathematical::pi*r/2));
        // } else {
        //     phi = 0;
        // }

        dirac *= phi;
    }

    // Return the final Dirac delta value
    dd = dirac/(h*h*h);
    return dd;
}
