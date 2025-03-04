use std::{fs::File, path::PathBuf};

use clap::{Parser, Subcommand};

use serde_json;

use petgraph::dot::{Dot, Config};

use pombase_gocam::{gocam_parse, make_gocam_model, GoCamRawModel};
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
    PrintActivities {
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
    },
    #[command(arg_required_else_help = true)]
    ConnectedGenes {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    NodesWithActivities {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    DetachedGenes {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    Serialize {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    }
}

fn print_tuples(model: &GoCamRawModel) {
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
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        model.id(),
        model.title(),
        subject.id,
        subject_type.label.as_ref().unwrap_or(empty),
        subject_type.id.as_ref().unwrap_or(&subject_type.type_string),
        fact.property_label,
        object.id,
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
                let model = make_gocam_model(&mut source)?;

                let stats = get_stats(&model);

                println!("{}\t{}\t{}\t{}\t{}\t{}\t{}", model.id(), model.taxon(),
                stats.total_genes, stats.total_complexes, stats.max_connected_activities,
                stats.total_connected_activities, stats.number_of_holes);
            }
        }
        Action::ConnectedGenes { paths } => {
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = make_gocam_model(&mut source)?;

                for (taxon, gene) in get_connected_genes(&model, 2) {
                    println!("{taxon}\t{gene}");
                }
            }
        }
        Action::NodesWithActivities { paths } => {
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = make_gocam_model(&mut source)?;

                for (taxon, gene) in get_connected_genes(&model, 1) {
                    println!("{taxon}\t{gene}");
                }
            }
        }
        Action::PrintTuples { paths } => {
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = gocam_parse(&mut source)?;
                print_tuples(&model);
            }
        },
        Action::PrintActivities { paths } => {
            println!("model_id\tmodel_title\ttaxon\tnode_id\tnode_label\tnode_type\tenabled_by_type\tenabled_by_id\tenabled_by_label\tprevious_nodes\tnext_nodes\tprocess\tinput\toutput\toccurs_in\tlocated_in");

            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = make_gocam_model(&mut source)?;

                let model_id = model.id();
                let model_title = model.title();
                let model_taxon = model.taxon();

                for node in model.node_iterator() {
                    println!("{}\t{}\t{}\t{}", model_id, model_title, model_taxon, node);
                }
            }
        },

        Action::FindHoles { paths } => {
            println!("model_id\tmodel_title\ttaxon\tactivity_id\tactivity_label\tprocess\tinput\toutput\toccurs_in\tlocated_in\ttype");
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = make_gocam_model(&mut source)?;

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

            let elements = model_to_cytoscape(&model);
            let elements_string = serde_json::to_string(&elements).unwrap();

            println!("{}", elements_string);
        }
        Action::CytoscapeSimple { path } => {
            let mut source = File::open(path).unwrap();
            let model = make_gocam_model(&mut source)?;

            let elements = model_to_cytoscape_simple(&model);
            let elements_string = serde_json::to_string(&elements).unwrap();

            println!("{}", elements_string);
        },
        Action::GraphVizDot { path } => {
            let mut source = File::open(path).unwrap();
            let model = make_gocam_model(&mut source)?;

            let dag_graphviz = Dot::with_attr_getters(
                model.graph(),
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
        },
        Action::DetachedGenes { paths } => {
//            println!("model_id\tmodel_title\ttaxon\tactivity_id\tactivity_label\tprocess\tinput\toutput\toccurs_in\tlocated_in\ttype");

            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = gocam_parse(&mut source)?;

                let model_id = model.id();

                let detached_genes = find_detached_genes(&model);

                for (id, gene_id, gene_label) in detached_genes {
                    println!("{}\t{}\t{}\t{}", model_id, id, gene_id, gene_label);

                }
            }
        },
        Action::Serialize { paths } => {
            let models: Vec<_> = paths.iter().map(|path| {
                let mut source = File::open(path).unwrap();
                let model = make_gocam_model(&mut source).unwrap();
                model
            })
            .collect();

            let models_string = serde_json::to_string(&models).unwrap();

            print!("{}", models_string);
        }
    }

    Ok(())
}
