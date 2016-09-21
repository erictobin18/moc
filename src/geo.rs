use std::cmp::Ordering;

pub trait HasLocation { // object exists at a single point (in the 2D axisymmetric plane), (r,z) : (f64,f64)
    fn location(&self) -> Point;
}

pub trait HasPath { // object exists as a path, a 2D countour (in the 2D axisymmetric plane), defined by a set of points it passes through
    fn point_list(&self) -> &Vec<Point>;
    fn r_at(&self, z: f64) -> f64; // interpolate the points to get the corresponding r value to the given z
    fn dr_dz_at(&self, z: f64) -> f64; // get the corresponding derivative (in z) of r at the given z
    fn intersect<T: HasPath>(&self, &other: T, guess: f64) -> Point; // find the intersection of two paths from an initial guess
    fn extrap_intersect<T: HasPath>(&self, self_z: f64, &other: T, other_z: f64) -> Point; // extrapolate two paths from initial points to find their intersection
    
}

#[derive(Copy,Clone,PartialEq)]
pub struct Point {
    pub r: f64, // distance from axis of symmetry
    pub z: f64, // distance from the most upstream point in the problem domain; +z is in flow direction
}

impl HasLocation for Point {
    fn location(&self) -> Point {
        *self
    }
}

impl Eq for Point {
    
}

impl PartialOrd for Point { 
    fn partial_cmp(&self, other: &Point) -> Option<Ordering> {
        if self.z < other.z {
            Some(Ordering::Less)
        }
        else if self.z > other.z {
            Some(Ordering::Greater)
        }
        else if self.r < other.r {
            Some(Ordering::Less)
        }
        else if self.r > other.r {
            Some(Ordering::Greater)
        }
        else {
            Some(Ordering::Equal)
        }
    }
}

impl Ord for Point {
    fn cmp(&self, other: &Point) -> Ordering {
        match self.partial_cmp(other) {
            Some(x) => x,
            None => Ordering::Equal,
        }
    }
}

#[derive(Clone)]
pub struct Path {
    // A path is a cubic spline interpolation between points
    pub points: Vec<Point>, // points the spline goes through
    pub dzi: Vec<f64>, // z[i + 1] - z[i]
    pub fpp: Vec<f64>, // second derivative at grid points
    pub generated: bool,     
}

#[allow(dead_code,unused_variables)]
impl Path {
    pub fn new(pts:Vec<Point>) -> Path {
        let n = pts.len();
        if n == 0 {
            Path { points: pts, dzi: vec![0.0;n], fpp: vec![0.0;n], generated: false }
        } else {
            Path { points: pts, dzi: vec![0.0;n-1], fpp: vec![0.0;n], generated: false }
        }
    }
    pub fn intersect(&self, other: &Path, guess: f64) -> Point {
        // search for root of other.r_at(z) - self.r_at(z)
        Point {z:0.0,r:0.0}
    }
    
    pub fn evaluate(&mut self, z:Point) -> f64 {
        if !self.generated {
            self.gen_spline();
        }
        assert!((self.points[0].z <= z.z) && (z.z <= self.points[self.points.len() - 1].z));
        let i = self.index(z);
        let dzp = self.points[i+1].z - z.z;
        let dzm = z.z - self.points[i].z;
        self.fpp[i]*(dzp*dzp*dzp/self.dzi[i] - self.dzi[i]*dzp)/6.0 +
            self.fpp[i+1]*(dzm*dzm*dzm/self.dzi[i] - self.dzi[i]*dzm)/6.0 +
            self.points[i].r*dzp/self.dzi[i] + 
            self.points[i+1].r*dzm/self.dzi[i]
    }
    pub fn evaluate_der(&mut self, z:Point) -> f64 {
        if !self.generated {
            self.gen_spline();
        }
        assert!((self.points[0].z <= z.z) && (z.z <= self.points[self.points.len() - 1].z));
        let i = self.index(z);
        let dzp = self.points[i+1].z - z.z;
        let dzm = z.z - self.points[i].z;

        0.0
    }
    // pub fn arc_len_point(&self, arc: f64) -> Point{
    //     if arc <= 0.0 { self.points[0] }
    //     else if arc >= 1.0 { self.points[self.points.len()-1] }
    //     else {
    //         let temp = match self.arc_lens.binary_search(&arc) {
    //             Ok(x)  => x,
    //             Err(x) => x,
    //         }
            
    //         // use the arc length formula for a cubic function to interpolate between temp and the arc length of the point on the other side of the interval


    //     }
    
    fn index(& self, z:Point) -> usize {
        let temp = match self.points.binary_search(&z) {
            Ok(x)  => x,
            Err(x) => x,
        };
        if temp <= 0 { 0 }
//        else if temp - 1 >= self.points.len() { self.points.len() - 1 }
        else { temp - 1 }
    }
    fn gen_spline(&mut self) {
        self.generated = true;
        // the piecewise expression to evaluate a spline depends on zi and ri, on dzi, and on fpp. To "generate" the spline, calculate the values of fpp from zi and ri
        // To find the values of fpp requires solving a system of N equations
        // The system is tridiagonal diagonal dominant, so it is unconditionally stable
        // a row of the matrix equations looks like:
        // {..., 0, 0, 0, dz_i-1/6, (dz_i-1 + dz_i)/3, dz_i/6, 0, 0, 0, ...}
        // the first row is {1, 0, 0, ...}, and the last is {..., 0, 0, 1}
        // the unknown vector is {fpp_0, fpp_1, ..., fpp_i, ..., fpp_n}
        // the constant vector is {0, (r_2 - r_1)/dz_1 - (r_1 - r_0)/dz_0, ..., (r_i+1 - r_i)/dz_i - (r_i - r_i-1)/dz_i-1, ..., 0}
        
        //self.points.sort_by(|a, b| a.partial_cmp(b).unwrap()); // floats don't have total ordering; must be sorted this way
        self.points.sort();
        self.dzi = self.points.iter().skip(1) // scan will take the difference of the ith and i-1th element; must skip the first
            .scan(self.points[0].z, |prev_z, pt| {
                let out = Some(pt.z - *prev_z); // dz[i] = z[i+1] - z[i], note dzi.len() == points.len() - 1
                *prev_z = pt.z;
                out
            } ).collect();
        // To solve a tridiagonal matrix:
        // [ b0 c0 0  0  0 ... ]
        // [ a1 b1 c1 0  0 ... ]
        // [ 0  a2 b2 c2 0 ... ]
        // [... 0  0  an-1 bn-1 cn-1 ]
        let n:usize = self.points.len();

        let a:Vec<f64> = (0..n).map(|x| {
            if x == 0 {
                0.0
            }
            else if x == n-1 {
                0.0
            }
            else {
                self.dzi[x-1]/6.0
            }
        } ).collect::<Vec<f64>>();

        let b:Vec<f64> = (0..n).map(|x| {
            if x == 0 {
                1.0
            }
            else if x == n-1 {
                1.0
            }
            else {
                (self.dzi[x-1] + self.dzi[x])/3.0
            }
        } ).collect::<Vec<f64>>();

        let c:Vec<f64> = (0..n).map(|x| {
            if x == 0 {
                0.0
            }
            else if x == n-1 {
                0.0
            }
            else {
                self.dzi[x]/6.0
            }
        } ).collect::<Vec<f64>>();

        let d:Vec<f64> = (0..n).map(|x| {
            if x == 0 {
                0.0
            }
            else if x == n-1 {
                0.0
            }
            else {
                (self.points[x+1].r - self.points[x].r)/self.dzi[x] - (self.points[x].r - self.points[x-1].r)/self.dzi[x-1]
            }
        } ).collect::<Vec<f64>>();

        let cp:Vec<f64> = c.iter().enumerate().scan(0.0,|prev_cp,x| {
            let out:f64 = x.1/(b[x.0] - a[x.0]*(*prev_cp as f64));
            *prev_cp = out;
            Some(out)
        } ).collect::<Vec<f64>>();

        let dp:Vec<f64> = d.iter().enumerate().scan(0.0,|prev_dp,x| {
            let out:f64 = (x.1-a[x.0]*(*prev_dp as f64))/(b[x.0] - a[x.0]*cp[x.0]);
            *prev_dp = out;
            Some(out)
        } ).collect::<Vec<f64>>();
        
        self.fpp = dp.iter().enumerate().rev().scan(0.0,|prev_fpp,x| {
            let out:f64 = x.1-cp[x.0]*(*prev_fpp as f64);
            *prev_fpp = out;
            Some(out)
        } ).collect::<Vec<f64>>();
    } 
}

impl Extend<Point> for Path {
    fn extend<T: IntoIterator<Item=Point>>(&mut self, iter:T) {
        self.points.extend(iter);
        self.generated = false;
    }

}

pub fn extrap_intersect(pt1: Point, pt1_der: f64, pt2: Point, pt2_der:f64)-> Point {
        // use first order appx for now
        let a = pt1_der;
        let b = pt2_der;

        let c = pt1.r - a*pt1.z;
        let d = pt2.r - b*pt2.z;
        
        Point {z:(d-c)/(a-b), r:(a*d-b*c)/(a-b)}
}

// solve y = a x + c, y = b x + d for x and y
pub fn solve_lin(a:f64, b:f64, c:f64, d:f64, x:&mut f64, y:&mut f64) {
    *x = (d-c)/(a-b);
    *y = (a*d - b*c)/(a-b);
}
