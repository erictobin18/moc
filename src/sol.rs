use geo;
use flu;

use std::fs::File;
use std::path::Path;
use std::io::Write;
use std::cmp;

#[derive(Clone)]
pub struct Solution {
    pub nodes: Vec<Node>,
    pub bounds: Vec<Boundary>,
    pub plus_char: Vec<CharLine>,
    pub minus_char: Vec<CharLine>,
    pub n: usize,
    pub stag_p: f64,
    pub stag_t: f64,
}

impl Solution {
    // pub fn new(bounds:Vec<Boundary>, stag_p:f64, stag_t:f64, num_init_nodes:usize) -> Solution {

    //     let mut node_vec:Vec<Node> = Vec::with_capacity(num_init_nodes);
    //     let mut pos_vec:Vec<CharLine> = Vec::with_capacity(num_init_nodes*2);
    //     let mut neg_vec:Vec<CharLine> = Vec::with_capacity(num_init_nodes);
    //     let mut pts:Vec<geo::Point> = Vec::new();
    //     let mut mach:f64 = 0.0;
    //     bounds.iter().find(|&&x| match x {
    //         Boundary::InitialLine{path:p,m:m,th:th} => {
    //             pts = p.points;
    //             mach = m[0];
    //             true
    //         },
    //         _ => false,
    //     } );
                
    //     let dr:f64 = pts[pts.len()-1].r/((num_init_nodes - 1) as f64);
    //     for i in 0..num_init_nodes {
    //         // for now, assume that the init line is vertical with a vertex at 0,0, the other vertex at the r-value of the last point in its path, and (m & th) constant and equal to the first value in their vectors
    //         node_vec.push( Node {
    //             pos: geo::Point { r: (i as f64)*dr, z: 0.0 },
    //             plus: num_init_nodes - i - 1,
    //             minus: i,
    //             state: flu::FluidState::default_state(mach,stag_p,stag_t),
    //             is_boundary_node: true,
                
    //         } );
    //         neg_vec.push( CharLine {
                
    //         } );
    //     }
    //     for i in 0..num_init_nodes {
    //         pos_vec.push( CharLine {
    //             is_plus_char: true,
    //             char_num: i,
    //             nodes: vec!(node_vec
                
    //         } );

    //     }
        
    //     Solution {
    //         nodes: node_vec,
    //         bounds: bounds,
    //         plus_char: pos_vec,
    //         minus_char: neg_vec,
    //         n: num_init_nodes - 1,
    //         stag_p: stag_p,
    //         stag_t: stag_t,
    //     }
    // }
    
    pub fn node_at(&self, plus:usize, minus:usize) -> &Node {
        if plus+self.n < minus || plus > minus+self.n || plus+minus < self.n { panic!() } // solution domain
        if minus + 1 > self.minus_char.len() { panic!() }

        &self.nodes[self.minus_char[minus].get_node_at_index(plus)]
    }
    pub fn print(&self) {
        let node_path = Path::new("sols/nodes.txt");
        let char_path = Path::new("sols/chars.txt");
        let info_path = Path::new("sols/info.txt");

        let mut node_file = match File::create(&node_path) {
            Err(why) => panic!(),
            Ok(file) => file,
        };
        let mut char_file = match File::create(&char_path) {
            Err(why) => panic!(),
            Ok(file) => file,
        };
        let mut info_file = match File::create(&info_path) {
            Err(why) => panic!(),
            Ok(file) => file,
        };
        write!(&mut node_file, "Node list:\nindex\tr\tz\tplus\tminus\tm\tth\n");
        for i in 0..self.nodes.len() {
            let nd = self.nodes[i].clone();
            write!(&mut node_file, "{}\t{:.*}\t{:.*}\t{}\t{}\t{:.*}\t{:.*}\n",i,4,nd.pos.r,4,nd.pos.z,nd.plus,nd.minus,4,nd.state.m,4,nd.state.th);
        }
        write!(&mut char_file, "Plus Char List:\nindex\tnodes\n");
        for i in 0..self.plus_char.len() {
            write!(&mut char_file, "{}\t",i);
            for j in 0..self.plus_char[i].nodes.len()-1 {
                write!(&mut char_file, "{},",self.plus_char[i].nodes[j]);
            }
            write!(&mut char_file, "{}\n",self.plus_char[i].nodes[self.plus_char[i].nodes.len() - 1]);
        }
        
        write!(&mut char_file, "Minus Char List:\nindex\tnodes\n");
        for i in 0..self.minus_char.len() {
            write!(&mut char_file, "{}\t",i);
            for j in 0..self.minus_char[i].nodes.len()-1 {
                write!(&mut char_file, "{},",self.minus_char[i].nodes[j]);
            }
            write!(&mut char_file, "{}\n",self.minus_char[i].nodes[self.minus_char[i].nodes.len() - 1]);
        }
        
        write!(&mut info_file, "Number of nodes: {}\n",self.nodes.len());
        write!(&mut info_file, "Number of init_line nodes: {}\n",self.n+1);
        write!(&mut info_file, "Stagnation pressure: {}\n",self.stag_p);
        write!(&mut info_file, "Stagnation temperature: {}\n",self.stag_t);
        write!(&mut info_file, "Inlet radius: {}\n",self.nodes[self.n].pos.r);
    }
    
    pub fn solve(&mut self) {
        let mut current_index:usize;        
        

        // initial line
        for j in 1..self.n+1 { // count from 1 to n-1, inclusive
            current_index = self.nodes.len();
            // j is the char_number of the - characteristic

            for i in (self.n-j+1)..(self.n+j) { // i is the + char #
                // create the node representing the intersection between the (+,-) characteristic (i,j)
                let new_node = flu::find_node(self.node_at(i,j-1),self.node_at(i-1,j));
                self.nodes.push(new_node);
                current_index = self.nodes.len();
                self.minus_char[j].nodes.push(current_index - 1);
                self.minus_char[j].path.extend(vec![self.nodes[current_index - 1].pos].into_iter());
                self.plus_char[i].nodes.push(current_index - 1);
                self.plus_char[i].path.extend(vec![self.nodes[current_index - 1].pos].into_iter());
            }
            // Centerline boundary node
            let new_node = flu::find_node_bnd_pls(self.node_at(self.n+j-1,j));
            self.nodes.push(new_node);
            current_index = self.nodes.len();
            self.minus_char[j].nodes.push(current_index - 1);
            self.minus_char[j].path.extend(vec![self.nodes[current_index - 1].pos].into_iter());
            self.plus_char[self.n + j].nodes.push(current_index - 1);
            self.plus_char[self.n + j].path.extend(vec![self.nodes[current_index - 1].pos].into_iter());

        }
    }
}

#[derive(Clone)]
pub struct Node {
    // container for all data associated with an intersection of a + and a - characteristic line
    pub pos: geo::Point,
    pub plus: usize,
    pub minus: usize,
    pub state: flu::FluidState,
    pub is_boundary_node: bool,
}

impl geo::HasLocation for Node {
    fn location(&self) -> geo::Point {
        self.pos
    }
}

#[allow(dead_code)]
#[derive(Clone)]
pub struct CharLine {
    // a + or a - characteristic line
    pub is_plus_char: bool,
    pub char_num: usize,
    pub nodes: Vec<usize>,
    pub integration_const: f64,
    pub path: geo::Path,
    pub char_index_offset: usize,
}

impl CharLine {
    fn get_node_at_index(&self, other_char:usize) -> usize {
        // get the node at intersection with other_char
        
        // first, determine the char number that goes through the first node in self                
        self.nodes[other_char - self.char_index_offset]        
    }
}

#[allow(dead_code)]
#[derive(Clone)]
pub enum Boundary {
    SymmetryAxis { point_a: geo::Point, point_b: geo::Point },
    Wall { path: geo::Path },
    InitialLine {
        path:geo::Path,
        m: Vec<f64>,
        th: Vec<f64>
    },
    Outlet { path: geo::Path },
    Error,
}
