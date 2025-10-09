use std::{collections::{BTreeSet, HashMap, HashSet}, fs::File, path::PathBuf};

use clap::{Parser, Subcommand};

use serde_json;

use petgraph::dot::{Dot, Config};

use pombase_gocam::{gocam_py::gocam_py_parse, parse_gocam_model,
                    raw::{gocam_parse_raw, GoCamRawModel}, GoCamEnabledBy, GoCamModel,
                    GoCamModelId, GoCamNode, GoCamNodeOverlap, GoCamNodeType, RemoveType};
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
    PrintNodes {
        #[arg(long)]
        remove_chemicals: bool,
        #[arg(long)]
        remove_inputs_outputs: bool,
        #[arg(required = true)]
        args: Vec<String>,
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
    CytoscapeSimpleMerged {
        #[arg(short, long)]
        taxon_id: Option<String>,
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    CytoscapeModelConnections {
        #[arg(short, long)]
        taxon_id: Option<String>,
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    CytoscapeModelConnectionsWithRelNodes {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
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
    AllGenes {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    GenesEnablingActivities {
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
    },
    #[command(arg_required_else_help = true)]
    OverlappingNodes {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    MakeChadoData {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    GocamPyParseTest {
        paths: Vec<PathBuf>,
    },
    #[command(arg_required_else_help = true)]
    JoiningChemicals {
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
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        model.id(),
        model.title(),
        model.taxon(),
        subject.id,
        subject_type.label.as_ref().unwrap_or(empty),
        subject_type.id.as_ref().unwrap_or(&subject_type.type_string),
        fact.property_label,
        object.id,
        object_type.label.as_ref().unwrap_or(empty),
        object_type.id.as_ref().unwrap_or(&object_type.type_string));
    }
}

fn node_type_summary_strings(node: &GoCamNode)
     -> (&str, &str, &str, String)
{
    match &node.node_type {
        GoCamNodeType::Unknown => ("unknown", "unknown", "unknown", "unknown".to_owned()),
        GoCamNodeType::Chemical => ("chemical", "", "", "".to_owned()),
        GoCamNodeType::UnknownMRNA => ("unknown_mrna", "", "", "".to_owned()),
        GoCamNodeType::MRNA(_) => ("mRNA", "", "", "".to_owned()),
        GoCamNodeType::Gene(_) => ("gene", "", "", "".to_owned()),
        GoCamNodeType::ModifiedProtein(_) => ("modified_protein", "", "", "".to_owned()),
        GoCamNodeType::Activity { enabler, .. } => match enabler {
            GoCamEnabledBy::Chemical(chem) => ("activity", "chemical", chem.id(), chem.label().to_owned()),
            GoCamEnabledBy::Gene(gene) => ("activity", "gene", gene.id(), gene.label()),
            GoCamEnabledBy::ModifiedProtein(prot) => ("activity", "modified_protein", prot.id(), prot.label().to_owned()),
            GoCamEnabledBy::Complex(complex) => ("activity", "complex", complex.id(), complex.label().to_owned()),
        }
    }
}


fn node_as_tsv(node: &GoCamNode) -> String {
    let mut ret = String::new();

    if let Some(ref original_model_id) = node.original_model_id {
        ret.push_str(original_model_id);
    }
    ret.push_str("\t");
    ret.push_str(&node.individual_gocam_id);
    ret.push_str("\t");

    ret.push_str(&format!("{}\t", node.node_id));

    ret.push_str(&format!("{}\t", node.label));

    let (node_type, enabled_by_type, enabled_by_id, enabled_by_label) =
        node_type_summary_strings(node);

    ret.push_str(&format!("{}\t{}\t{}\t{}\t", node_type, enabled_by_type, enabled_by_id, enabled_by_label));

    if let Some(ref part_of_process) = node.part_of_process {
        ret.push_str(&format!("{}\t", part_of_process.label_or_id()));
    } else {
        ret.push_str(&format!("\t"));
    }

    if let GoCamNodeType::Activity { ref inputs, ref outputs, .. } = node.node_type {
    let has_input_string =
        inputs.iter().map(|l| l.to_string()).collect::<Vec<_>>().join(",");
    if has_input_string.len() > 0 {
        ret.push_str(&format!("{}\t", has_input_string));
    } else {
        ret.push_str(&format!("\t"));
    }
    let has_output_string =
        outputs.iter().map(|l| l.to_string()).collect::<Vec<_>>().join(",");
    if has_output_string.len() > 0 {
        ret.push_str(&format!("{}\t", has_output_string));
    } else {
        ret.push_str(&format!("\t"));
    }
} else {
    ret.push_str(&format!("\t"));
}

    let occurs_in_string = node.occurs_in
        .iter()
        .map(|occurs_in| occurs_in.id())
        .collect::<Vec<_>>()
        .join(",");
    ret.push_str(&format!("{}\t", occurs_in_string));

    if let Some(ref located_in) = node.located_in {
        ret.push_str(&format!("{}", located_in.label_or_id()));
    }
    ret.push_str(&format!("\t"));

    if let Some(ref happens_during) = node.happens_during {
        ret.push_str(&format!("{}", happens_during.label_or_id()));
    }
    ret.push_str(&format!("\t"));

    if let GoCamNodeType::Activity { ref enabler, .. } = node.node_type {
        if let GoCamEnabledBy::Complex(ref complex) = enabler {
            let parts = complex.has_part_genes.iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>().join(",");
            ret.push_str(&parts);
        }
    }

    ret
}

fn has_connected_genes(model: &GoCamModel) -> bool {
    let connected_genes_by_activity_count = get_connected_genes(&model);
    connected_genes_by_activity_count.get(&2).is_some()
}

fn filter_models_by_org(models: &[GoCamModel], taxon: &str)
    -> Vec<GoCamModel>
{
    return models.iter()
        .filter(|model| model.taxon().contains(taxon))
        .cloned()
        .collect()
}

fn models_from_paths(paths: &Vec<PathBuf>)
    -> Vec<GoCamModel>
{
    let models: Vec<_> = paths.iter().map(|path| {
        let mut source = File::open(path).unwrap();
        let model = parse_gocam_model(&mut source).unwrap();
        model
    }).collect();

    models.into_iter().collect()
}

fn model_from_paths(paths_string: &str)
    -> GoCamModel
{
    let paths: Vec<PathBuf> = paths_string.split('+').map(PathBuf::from).collect();
    let models = models_from_paths(&paths);

    if models.len() > 1 {
        GoCamModel::merge_models("merged", "merged models", &models).unwrap()
    } else {
        models.into_iter().next().unwrap()
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    match args.action {
        Action::Stats { paths } => {
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = parse_gocam_model(&mut source)?;

                let stats = get_stats(&model);

                println!("{}\t{}\t{}\t{}\t{}\t{}\t{}", model.id(), model.taxon(),
                stats.total_genes, stats.total_complexes, stats.max_connected_activities,
                stats.total_connected_activities, stats.number_of_holes);
            }
        }
        Action::ConnectedGenes { paths } => {
            println!("taxon\tgene");
            let mut seen_genes = HashSet::new();

            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = parse_gocam_model(&mut source)?;

                let connected_genes_by_activity_count =
                    get_connected_genes(&model);
                if let Some(connected_genes) = connected_genes_by_activity_count.get(&2) {
                    for gene in connected_genes {
                        if seen_genes.contains(&(model.taxon().to_owned(), gene.to_owned())) {
                            continue;
                        } else {
                            seen_genes.insert((model.taxon().to_owned(), gene.to_owned()));
                        }
                        println!("{}\t{gene}", model.taxon());
                    }
                }
            }
        }
        Action::AllGenes { paths } => {
            println!("taxon\tgene");
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = parse_gocam_model(&mut source)?;

                let connected_genes_by_activity_count =
                    get_connected_genes(&model);
                if let Some(connected_genes) = connected_genes_by_activity_count.get(&1) {
                    for gene in connected_genes {
                        println!("{}\t{gene}", model.taxon());
                    }
                }
            }
        }
        Action::GenesEnablingActivities { paths } => {
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = parse_gocam_model(&mut source)?;

                let genes = model.genes_enabling_activities();

                for (gene_id, _) in &genes {
                    println!("{gene_id}");
                }
            }
        }
        Action::PrintTuples { paths } => {
            println!("model_id\tmodel_title\ttaxon\tsubject_id\tsubject_label\tsubject_type_id\trelation\tobject_id\tobject_label\tobject_type_id");

            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = gocam_parse_raw(&mut source)?;
                print_tuples(&model);
            }
        },
        Action::PrintNodes { remove_chemicals, remove_inputs_outputs, args } => {
            println!("model_id\tmodel_title\ttaxon\toriginal_model_id\tindividual_gocam_id\tnode_id\tnode_label\tnode_type\tenabled_by_type\tenabled_by_id\tenabled_by_label\tprocess\tinput\toutput\toccurs_in\tlocated_in\thappens_during\tparts");

            for arg in args {

                let model = {
                    let model = model_from_paths(&arg);

                    let mut remove_types = HashSet::new();

                    if remove_chemicals {
                        remove_types.insert(RemoveType::Chemicals);
                    } else {
                        if remove_inputs_outputs {
                            remove_types.insert(RemoveType::Targets);
                        }
                    }

                    if remove_types.is_empty() {
                        model
                    } else {
                        model.remove_nodes(remove_types)
                    }
                };

                let model_id = model.id();
                let model_title = model.title();
                let model_taxon = model.taxon();

                for (_, node) in model.node_iterator() {
                    println!("{}\t{}\t{}\t{}", model_id, model_title, model_taxon, node_as_tsv(node));
                }
            }
        },

        Action::FindHoles { paths } => {
            println!("model_id\tmodel_title\ttaxon\toriginal_model_id\tindividual_gocam_id\tnode_id\tnode_label\tnode_type\tenabled_by_type\tenabled_by_id\tenabled_by_label\tprocess\tinput\toutput\toccurs_in\tlocated_in\thappens_during\tparts");
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = parse_gocam_model(&mut source)?;

                let model_id = model.id();
                let model_title = model.title();
                let model_taxon = model.taxon();

                let hole_nodes = find_holes(&model);

                for hole_node in hole_nodes {
                    println!("{}\t{}\t{}\t{}", model_id, model_title, model_taxon,
                             node_as_tsv(&hole_node));

                }
            }
        },
        Action::Cytoscape { path } => {
            let mut source = File::open(path).unwrap();
            let model = gocam_parse_raw(&mut source)?;

            let elements = model_to_cytoscape(&model);
            let elements_string = serde_json::to_string(&elements).unwrap();

            println!("{}", elements_string);
        }
        Action::CytoscapeSimple { path } => {
            let mut source = File::open(path).unwrap();
            let model = parse_gocam_model(&mut source)?;

            let elements = model_to_cytoscape_simple(&model, &vec![],
                                                     GoCamCytoscapeStyle::IncludeParents);
            let elements_string = serde_json::to_string(&elements).unwrap();

            println!("{}", elements_string);
        },
        Action::CytoscapeSimpleMerged { taxon_id, paths } => {
            let models: Vec<_> =
                if let Some(taxon_id) = taxon_id {
                    let taxon_id = taxon_id.strip_prefix("NCBITaxon:").unwrap_or(&taxon_id);
                    filter_models_by_org(&models_from_paths(&paths), taxon_id)
                } else {
                    models_from_paths(&paths)
                }
                .into_iter().filter(has_connected_genes).collect();
            let merged = GoCamModel::merge_models("merged", "merged models", &models)?;

            let elements = model_to_cytoscape_simple(&merged, &vec![], GoCamCytoscapeStyle::IncludeParents);
            let elements_string = serde_json::to_string(&elements).unwrap();

            println!("{}", elements_string);
        },
        Action::CytoscapeModelConnections { taxon_id, paths } => {
            let all_models = models_from_paths(&paths);
            let models: Vec<_> =
                if let Some(taxon_id) = taxon_id {
                    let taxon_id = taxon_id.strip_prefix("NCBITaxon:").unwrap_or(&taxon_id);
                    filter_models_by_org(&all_models, taxon_id)
                } else {
                    models_from_paths(&paths)
                }
                .into_iter().filter(has_connected_genes).collect();

            let overlaps = GoCamModel::find_overlaps(&models);

            let model_ids_and_titles: Vec<_> =
                all_models.iter()
                .map(|model| (model.id().to_owned(), model.title().to_owned()))
                .collect();
            let elements = model_connections_to_cytoscope(&overlaps, &model_ids_and_titles);

            let elements_string = serde_json::to_string(&elements).unwrap();

            println!("{}", elements_string);
        },
        Action::CytoscapeModelConnectionsWithRelNodes { paths } => {
            let models = models_from_paths(&paths);

            let elements = model_pathways_to_cytoscope_test(&models);

            let elements_string = serde_json::to_string(&elements).unwrap();

            println!("{}", elements_string);
        },
        Action::GraphVizDot { path } => {
            let mut source = File::open(path).unwrap();
            let model = parse_gocam_model(&mut source)?;

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
                let model = gocam_parse_raw(&mut source)?;

                let model_id = model.id();

                let detached_genes = find_detached_genes(&model);

                for (id, gene_id, gene_label) in detached_genes {
                    println!("{}\t{}\t{}\t{}", model_id, id, gene_id, gene_label);

                }
            }
        },
        Action::Serialize { paths } => {
            let models = models_from_paths(&paths);

            let models_string = serde_json::to_string(&models).unwrap();

            print!("{}", models_string);
        },
        Action::OverlappingNodes { paths } => {
            let models = models_from_paths(&paths);

            let overlaps = GoCamModel::find_overlaps(&models);

            println!("model_titles\tmodel_ids\tid\tlabel\tdescription\tpart_of_process\toccurs_in\tlocated_in");

            for overlap in &overlaps {
                let mut models = HashSet::new();
                models.extend(&overlap.models);
                if models.len() == 1 {
                    continue;
                }

                let process_label =
                    if let Some(ref part_of_process) = overlap.part_of_process {
                        part_of_process.label.clone()
                    } else {
                        String::default()
                    };

                let occurs_in_label = overlap.occurs_in
                    .iter()
                    .map(|occurs_in| occurs_in.label())
                    .collect::<Vec<_>>()
                    .join(",");

                let located_in_label =
                    if let Some(ref located_in) = overlap.located_in {
                        located_in.label().to_owned()
                    } else {
                        String::default()
                    };

                let (model_ids, model_titles): (Vec<_>, Vec<_>) =
                    overlap.models.iter()
                           .map(|(id, title, _)| (id.to_owned(), title.to_owned()))
                           .unzip();

                    println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                         model_titles.into_iter().collect::<Vec<_>>().join(","),
                         model_ids.into_iter().collect::<Vec<_>>().join("+"),
                         overlap.node_id, overlap.node_label,
                         overlap.node_type,
                         process_label,
                         occurs_in_label,
                         located_in_label);
            }
        },
        Action::MakeChadoData { paths } => {
            let models = models_from_paths(&paths);

            let data_for_chado = make_chado_data(&models);

            let chado_string = serde_json::to_string(&data_for_chado).unwrap();

            println!("{}", chado_string);
        },
        Action::GocamPyParseTest { paths } => {
            for path in paths {
                let mut source = File::open(path).unwrap();
                let gocam_py_model = gocam_py_parse(&mut source)?;
                println!("id: {}", gocam_py_model.id);
            }
        },
        Action::JoiningChemicals { paths } => {
            let models = models_from_paths(&paths);

            let overlaps = GoCamModel::find_overlaps(&models);

            let overlaps_by_models: HashMap<BTreeSet<GoCamModelId>, GoCamNodeOverlap> = overlaps
                .into_iter()
                .map(|overlap| {
                    let ids = overlap.models.iter().map(|(id,_, _)| id.to_owned()).collect();
                    (ids, overlap)
                })
                .collect();

            let mut chemical_groups = HashMap::new();

            for model in models {
                for (_, node) in model.node_iterator() {
                    if node.node_type == GoCamNodeType::Chemical {
                        println!("{}", node);
                        let (model_id, _) = node.models.first().unwrap();
                        chemical_groups.entry((node.node_id.clone(), node.label.clone()))
                            .or_insert_with(BTreeSet::new)
                            .insert(model_id.to_owned());
                    }
                }
            }

            for ((chem_id, chem_label), model_ids) in chemical_groups.into_iter() {
                if model_ids.len() < 2 {
                    continue;
                }

                if overlaps_by_models.contains_key(&model_ids) {
                    continue;
                }

                let model_ids_string = model_ids.into_iter().collect::<Vec<_>>().join(",");

                println!("{}\t{}\t{}", chem_id, chem_label, model_ids_string);
            }
        }
    }

    Ok(())
}
