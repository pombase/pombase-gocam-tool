use std::{collections::{HashMap, HashSet}, fmt::{self, Display}, fs::File, path::PathBuf};

use clap::{Parser, Subcommand};
extern crate serde_json;
#[macro_use] extern crate serde_derive;

use petgraph::{dot::Config, graph::NodeIndex, visit::{EdgeRef, IntoNodeReferences, NodeRef}, Graph, Undirected};
use petgraph::visit::Bfs;

use pombase_gocam::{gocam_parse, FactId, GoCamModel, Individual, IndividualId, IndividualType};
use pombase_gocam_process::*;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Args {
    #[command(subcommand)]
    action: Action,
}

#[derive(Subcommand)]
enum Action {
    #[command(arg_required_else_help = true)]
    Stats {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    PrintTuples {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    GraphTest {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    FindHoles {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    Cytoscape {
        #[arg(required = true)]
        path: PathBuf,
    },
    #[command(arg_required_else_help = true)]
    CytoscapeSimple {
        #[arg(required = true)]
        path: PathBuf,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    match args.action {
        Action::Stats { paths } => {
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = gocam_parse(&mut source)?;

                print_stats(&model);
            }
        }
        Action::PrintTuples { paths } => {
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = gocam_parse(&mut source)?;
                print_tuples(&model);
            }
        },
        Action::GraphTest { paths } => {
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = gocam_parse(&mut source)?;
                let graph = make_graph(&model);

                let first_index = graph.node_indices().next().unwrap();

                let mut bfs = Bfs::new(&graph, first_index);
                while let Some(_) = bfs.next(&graph) {
                   //  println!("{}", model.id());
                }
            }
        },
        Action::FindHoles { paths } => {
            println!("model_id\tmodel_title\ttaxon\tactivity_id\tactivity_label\ttype\tprocess\tinput\toutput\toccurs_in\tlocated_in");
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = gocam_parse(&mut source)?;

                let model_id = model.id();
                let model_title = model.title();
                let model_taxon = model.taxon();

                let hole_nodes = find_holes(&model);

                for hole_node in hole_nodes {
                    println!("{}\t{}\t{}\t{}", model_id, model_title, model_taxon,
                             hole_node);

                }
            }
        },
        Action::Cytoscape { path } => {
            let mut source = File::open(path).unwrap();
            let model = gocam_parse(&mut source)?;

            let cytoscape_text = model_to_cytoscape(&model);

            println!("{}", cytoscape_text);
        }
        Action::CytoscapeSimple { path } => {
            let mut source = File::open(path).unwrap();
            let model = gocam_parse(&mut source)?;
            let graph = make_graph(&model);

            let cytoscape_text = model_to_cytoscape_simple(&graph);

            println!("{}", cytoscape_text);
        }
    }

    Ok(())
}
