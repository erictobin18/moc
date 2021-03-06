mod geo;
mod sol;
mod flu;


#[allow(unused_mut)]
fn main() {
    let mut nds:Vec<sol::Node> = Vec::new();

    let n:usize = 10;
    let dr:f64 = 0.001;
    let p0 = 5000000.0;
    let t0 = 4000.0;
    let mach = 2.0;

    let wall_path = geo::Path::new(vec![
        geo::Point { r: (n as f64)*dr, z: 0.0 },
        geo::Point { r: (n as f64)*dr, z: 100.0 },
    ]);
    let out_path = geo::Path::new(vec![
        geo::Point { r: (n as f64)*dr, z: 100.0 },
        geo::Point { r: 0.0, z: 100.0 },
    ]);

    let mut i_point:Vec<geo::Point> = Vec::new();
    let mut i_m:Vec<f64> = Vec::new();
    let mut i_th:Vec<f64> = Vec::new();
        
    for i in 0..n+1 {
        let pt = geo::Point { r: (i as f64)*dr, z: 0.0};
        nds.push( sol::Node {
            pos: pt,
            plus: n-i,
            minus: i,
            state: flu::FluidState::default_state(mach-(i as f64)*0.05,p0,t0),
            is_boundary_node: true,
        } );
        i_point.push(pt);
        i_m.push(mach);
        i_th.push(0.0);
    }

    let wall = sol::Boundary::Wall { path: wall_path };
    let out = sol::Boundary::Outlet { path: out_path };
    let axis = sol::Boundary::SymmetryAxis {
        point_a: geo::Point { r: 0.0, z: 0.0 },
        point_b: geo::Point { r: 0.0, z: 100.0 },
    };
    let i_line = sol::Boundary::InitialLine {
        path: geo::Path::new(i_point),
        m: i_m,
        th: i_th,
    };

    let bds:Vec<sol::Boundary> = vec![i_line, wall, out, axis];

    let mut pls_vec:Vec<sol::CharLine> = Vec::with_capacity(2*n + 1);
    let mut neg_vec:Vec<sol::CharLine> = Vec::with_capacity(n + 1);
    
    for i in 0..n+1 {
        pls_vec.push(sol::CharLine {
            is_plus_char: true,
            char_num: i,
            nodes: vec![n-i],
            integration_const: 0.0,
            path: geo::Path::new(vec![nds[n - i].pos]),
            char_index_offset: n - i,
        });
        neg_vec.push(sol::CharLine {
            is_plus_char: false,
            char_num: i,
            nodes: vec![i],
            integration_const: 0.0,
            path: geo::Path::new(vec![nds[i].pos]),
            char_index_offset: n - i,            
        });
    }
    for i in n+1..2*n+1 {
        pls_vec.push(sol::CharLine {
            is_plus_char: true,
            char_num: i,
            nodes: vec![],
            integration_const: 0.0,
            path: geo::Path::new(vec![]),
            char_index_offset: i - n,
        });
    }
    
    
    // Create a problem to give to the solver
    let mut problem = sol::Solution { nodes:nds, bounds:bds, plus_char: pls_vec, minus_char: neg_vec, n:n, stag_p: p0, stag_t: t0};
    
    problem.solve();

    problem.print();
}
