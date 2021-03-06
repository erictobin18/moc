use geo;
use sol;
//use std::f64;

const GAMMA:f64 = 1.4;

#[allow(dead_code, unused_variables)]
#[derive(Clone)]
pub struct FluidState {
    pub m: f64, // local Mach number; m = v / a, a = sqrt(GAMMA*specific_gas_const*temp)
    pub th: f64, // local flow angle; angle between v and +z axis; also known as turning angle
    vz: f64, // local axial velocity
    vr: f64, // local radial velocity
    temp: f64, // local temperature (static temperature is constant)
    p: f64, // local pressure (static pressure is constant)
    rho: f64, // local density (derived quantity from the equation of state, the ideal gas law)
}

impl FluidState {
    pub fn default_state(mach:f64, p0:f64, t0:f64) -> FluidState {
        let t = 0.0; // WRONG -- calculate from M and T0
        FluidState {
            m: mach,
            th: 0.0,
            vz: 0.0, // WRONG -- calculate from M and T
            vr: 0.0,
            temp: t,
            p: 0.0, // WRONG -- calculate from M and p0
            rho: 0.0 // WRONG -- calculate from ideal gas law
        }
    }
}

pub fn find_node(pos:&sol::Node, neg:&sol::Node) -> sol::Node {
    let pos_mu = f64::asin(1.0/pos.state.m);
    let neg_mu = f64::asin(1.0/neg.state.m);

    // use path data from nodes to find intersection
    let intersection = geo::extrap_intersect(pos.pos, f64::tan(pos.state.th + pos_mu), neg.pos, f64::tan(neg.state.th - neg_mu));
    
    
    let pos_dr = intersection.r - pos.pos.r;
    let neg_dr = intersection.r - neg.pos.r;
    
    // use fluid data from nodes to find fluid info
    // start the iterative analysis with averages
    let mut mach = (pos.state.m + neg.state.m)/2.0;
    let mut theta = (pos.state.th + neg.state.th)/2.0;

//    println!("intersection: r={}, z={}",intersection.r,intersection.z);
//    println!("m = {}, th = {}",mach,theta);
    let mut loop_index = 0;
    loop {
        loop_index = loop_index + 1;
        // process to solve simultaneous equations:
        // th - th1 - (1 - m1/m) sqrt(m^2-1)/(1+m*(GAMMA -1)/2) - dr1/r (1/(cot(th)-sqrt(m^2-1) == 0
        // th - th2 + (1 - m1/m) sqrt(m^2-1)/(1+m*(GAMMA -1)/2) - dr1/r (1/(cot(th)+sqrt(m^2-1) == 0
        // note the second equation has two different signs in it
        // *1 refers to the +char, *2 refers to the -char
        // define f(m, th) = eq1, g(m, th) = eq2
        // calculate partial derivatives of f & g w.r.t. m & th
        // set the { first-order Taylor expansion of f and g in delta_m and delta_th around m and th } equal to 0
        // solve for delta_m and delta_th
        // these values are appx. for the difference between the guessed values of m and th and their actual values
        // if delta_m and delta_th are small, return m+delta_m and th+delta_th as the correct values of m and th
        // otherwise, m = m + delta_m and th = th + delta_th are the new guesses, loop
        
        // calculate function values

//        println!("pos.state.th: {}",pos.state.th);
        let f = theta - pos.state.th
            - (pos_dr/intersection.r)/((1.0/f64::tan(theta)) - f64::sqrt(mach.powi(2) - 1.0))
            - (1.0 - pos.state.m/mach)*f64::sqrt(mach.powi(2) - 1.0)/(1.0 + (GAMMA - 1.0)*mach.powi(2)/2.0);
        let g = theta - neg.state.th
            - (neg_dr/intersection.r)/((1.0/f64::tan(theta)) + f64::sqrt(mach.powi(2) - 1.0))
            + (1.0 - neg.state.m/mach)*f64::sqrt(mach.powi(2) - 1.0)/(1.0 + (GAMMA - 1.0)*mach.powi(2)/2.0);

//        println!("f = {}, g = {}",f,g);
            
        // calculate partial derivatives
        let dfdm = -(pos_dr/intersection.r)*(mach/f64::sqrt(mach.powi(2) - 1.0))/(1.0/f64::tan(theta) - f64::sqrt(mach.powi(2) - 1.0)).powi(2)
            + (2.0*f64::sqrt(mach.powi(2) - 1.0)*(mach-pos.state.m)/(mach.powi(2) + 2.0/(GAMMA - 1.0)) - (mach + pos.state.m/mach.powi(2))/f64::sqrt(mach.powi(2) - 1.0))/(1.0 + (mach.powi(2))*(GAMMA - 1.0)/2.0);
        let dgdm = (neg_dr/intersection.r)*(mach/f64::sqrt(mach.powi(2) - 1.0))/(1.0/f64::tan(theta) + f64::sqrt(mach.powi(2) - 1.0)).powi(2)
            - (2.0*f64::sqrt(mach.powi(2) - 1.0)*(mach-neg.state.m)/(mach.powi(2) + 2.0/(GAMMA - 1.0)) - (mach + neg.state.m/mach.powi(2))/f64::sqrt(mach.powi(2) - 1.0))/(1.0 + (mach.powi(2))*(GAMMA - 1.0)/2.0);
        let dfdth = 1.0 - ((pos_dr/intersection.r)/(f64::cos(theta) - f64::sqrt(mach.powi(2) - 1.0))).powi(2);
        let dgdth = 1.0 - ((neg_dr/intersection.r)/(f64::cos(theta) + f64::sqrt(mach.powi(2) - 1.0))).powi(2);

//        println!("dfdm = {}, dgdm = {}, dfdth = {}, dgdth = {}",dfdm,dgdm,dfdth,dgdth);

        // delta_m == (-dfdth/dfdm)*delta_th-f(m,th)
        // delta_m == (-dgdth/dgdm)*delta_th-g(m,th)
        // a = -dfdth/dfdm, b = -dgdth/dgdm, c = -f(m,th), d = -g(m,th)
        let mut delta_m:f64 = 0.0;
        let mut delta_th:f64 = 0.0;

        geo::solve_lin(-dfdth/dfdm, -dgdth/dgdm, -f, -g, &mut delta_th, &mut delta_m);
//        println!("delta_m = {}, delta_th = {}",delta_m,delta_th);

//        println!("error m: {}", f64::abs(delta_m/mach));
//        println!("error th: {}", f64::abs(delta_th/theta));
        if f64::abs(delta_m/mach) < 1e-9 && (f64::abs(delta_th/theta) < 1e-9 || (theta == 0.0 && f64::abs(delta_th) < 1e-9)) {
            mach = mach + delta_m;
            theta = theta + delta_th;
            break;
        }
        mach = mach + delta_m;
        theta = theta + delta_th;

//        println!("m = {}, th = {}",mach,theta);
        if loop_index > 100 {
            println!("Find node solver loop iterated more than 100 times! Panic!ing");
            panic!();
        }                                   
    }
    
    sol::Node {
        pos: intersection, // unknown, 2 vars
        plus: pos.plus,
        minus: neg.minus,
        state: FluidState { // unknown, 2 vars
            m: mach,
            th: theta,
            vz: 0.0, // WRONG
            vr: 0.0, // WRONG
            temp: 0.0, // WRONG
            p: 0.0, // WRONG
            rho: 0.0 }, // WRONG               
        is_boundary_node: false,
        // only true if it's the last node, the one on the centerline
    }
}

pub fn find_node_bnd_pls(neg:&sol::Node) -> sol::Node {
     sol::Node {
        pos: geo::Point { r: 0.0, z: neg.pos.z + 0.8660254037844385 },
        plus: neg.plus+1,
        minus: neg.minus,
        state: FluidState {
            m: 2.0,
            th: 0.0,
            vz: 0.0,
            vr: 0.0,
            temp: 0.0,
            p: 0.0,
            rho: 0.0 },                    
        is_boundary_node: true,
        // only true if it's the last node, the one on the centerline
    }
}

#[allow(dead_code)]
pub fn find_node_bnd_neg(pos:&sol::Node) -> sol::Node {
     sol::Node {
        pos: geo::Point { r: 0.0, z: 0.0 },
        plus: pos.plus,
        minus: pos.minus + 1,
        state: FluidState {
            m: 0.0,
            th: 0.0,
            vz: 0.0,
            vr: 0.0,
            temp: 0.0,
            p: 0.0,
            rho: 0.0 },                    
        is_boundary_node: true,
        // only true if it's the last node, the one on the centerline
    }
}

