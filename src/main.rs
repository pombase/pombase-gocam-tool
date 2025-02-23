use std::{fs::File, path::PathBuf};

use clap::{Parser, Subcommand};

use petgraph::dot::{Dot, Config};
use petgraph::visit::Bfs;

use pombase_gocam::{gocam_parse, GoCamModel};
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
    #[command(arg_required_else_help = true)]
    GraphVizDot {
        #[arg(required = true)]
        path: PathBuf,
    }
}

fn print_tuples(model: &GoCamModel) {
    let empty = &"".to_owned();
    for fact in model.facts() {
        let subject = model.fact_subject(fact);
        let object = model.fact_object(fact);
        let Some(subject_type) = subject.types.get(0)
        else {
            continue;
        };
        let Some(object_type) = object.types.get(0)
        else {
            continue;
        };
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}",
                 model.id(),
                 model.title(),
                 subject_type.label.as_ref().unwrap_or(empty),
                 subject_type.id.as_ref().unwrap_or(&subject_type.type_string),
                 fact.property_label,
                 object_type.label.as_ref().unwrap_or(empty),
                 object_type.id.as_ref().unwrap_or(&object_type.type_string));
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    match args.action {
        Action::Stats { paths } => {
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = gocam_parse(&mut source)?;

                let stats = get_stats(&model);

                println!("{}\t{}\t{}\t{}\t{}\t{}", model.id(), model.taxon(),
                         stats.total_genes, stats.max_connected_genes,
                         stats.total_connected_genes, stats.number_of_holes);
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
                while let Some(nx) = bfs.next(&graph) {
                    let node = graph.node_weight(nx).unwrap();
                    println!("{}: {} {}", model.id(), node.id, node.label);
                }
            }
        },
        Action::FindHoles { paths } => {
            println!("model_id\tmodel_title\ttaxon\tactivity_id\tactivity_label\tprocess\tinput\toutput\toccurs_in\tlocated_in\ttype");
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
        },
        Action::GraphVizDot { path } => {
            let mut source = File::open(path).unwrap();
            let model = gocam_parse(&mut source)?;
            let graph = make_graph(&model);

            let dag_graphviz = Dot::with_attr_getters(
                &graph,
                &[Config::NodeNoLabel, Config::EdgeNoLabel],
                &|_, edge| format!("label = \"{}\"", edge.weight().label),
                &|_, (_, node)| {
                    let enabler_label = node.enabler_label();
                    if enabler_label.len() > 0 {
                        format!("label = \"{}\"", enabler_label)
                    } else {
                        format!("label = \"{}\"", node.label)
                    }
                },
            );

            println!("{}", dag_graphviz);
        }
    }

    Ok(())
}
