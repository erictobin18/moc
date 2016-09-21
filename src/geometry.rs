trait HasLocation { // object exists at a single point (in the 2D axisymmetric plane), (r,z) : (f64,f64)
    fn location(&self) -> Point;
}

trait HasPath { // object exists as a path, a 2D countour (in the 2D axisymmetric plane), defined by a set of points it passes through
    fn point_list(&self) -> Vec<Point>;
    fn r_at(&self, z: f64) -> f64; // interpolate the points to get the corresponding r value to the given z
    fn dr_dz_at(&self, z: f64) -> f64; // get the corresponding derivative (in z) of r at the given z
    fn intersect<T: HasPath>(&self, guess: Point) -> Point; // find the intersection of two paths from an initial guess
}

#[derive(Copy,Clone)]
struct Point {
    r: f64, // distance from axis of symmetry
    z: f64, // distance from the most upstream point in the problem domain; +z is in flow direction
}

impl HasLocation for Point {
    fn location(&self) -> Point {
        *self
    }
}

struct Path {
    // A path is a cubic spline interpolation between points
    points: Vec<Point>, // points the spline goes through
    dzi: Vec<f64>, // z[i + 1] - z[i]
    fpp: Vec<f64>, // second derivative at grid points
    
}

impl Path {
    fn gen_spline(&mut self) {
        // the piecewise expression to evaluate a spline depends on zi and ri, on dzi, and on fpp. To "generate" the spline, calculate the values of fpp from zi and ri
        // To find the values of fpp requires solving a system of N equations
        // The system is tridiagonal diagonal dominant, so it is unconditionally stable
        // a row of the matrix equations looks like:
        // {..., 0, 0, 0, dz_i-1/6, (dz_i-1 + dz_i)/3, dz_i/6, 0, 0, 0, ...}
        // the first row is {1, 0, 0, ...}, and the last is {..., 0, 0, 1}
        // the unknown vector is {fpp_0, fpp_1, ..., fpp_i, ..., fpp_n}
        // the constant vector is {0, (r_2 - r_1)/dz_1 - (r_1 - r_0)/dz_0, ..., (r_i+1 - r_i)/dz_i - (r_i - r_i-1)/dz_i-1, ..., 0}
        

        
    } 
}
