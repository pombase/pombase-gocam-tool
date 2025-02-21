use std::{collections::HashMap, fmt::{self, Display}, fs::File, path::PathBuf};

use clap::{Parser, Subcommand};

use petgraph::Graph;
use petgraph::visit::Bfs;

use pombase_gocam::{gocam_parse, FactId, GoCamModel, Individual, IndividualId, IndividualType};

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
}

type GoCamGraph = Graph::<GoCamNode, GoCamEdge>;

struct GoCamGraphModel {
    pub title: String,
    pub graph: GoCamGraph,
}

type GoCamComplex = IndividualType;
type GoCamGene = IndividualType;
type GoCamChemical = IndividualType;
type GoCamModifiedProtein = IndividualType;
type GoCamComponent = IndividualType;
type GoCamProcess = IndividualType;
type GoCamInput = IndividualType;
type GoCamOutput = IndividualType;

#[derive(Clone, Debug)]
enum GoCamActivity {
    Complex(GoCamComplex),
    Gene(GoCamGene),
    Chemical(GoCamChemical),
    ModifiedProtein(GoCamModifiedProtein),
}

impl GoCamActivity {
    pub fn label(&self) -> &str {
        let maybe_label = match self {
            GoCamActivity::Complex(complex) => &complex.label,
            GoCamActivity::Gene(gene) => &gene.label,
            GoCamActivity::Chemical(chemical) => &chemical.label,
            GoCamActivity::ModifiedProtein(modified_protein) => &modified_protein.label,
        };
        maybe_label.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN")
    }
}

#[derive(Clone, Debug)]
enum GoCamNodeType {
    Unknown,
    Chemical,
    Activity(GoCamActivity),
}

#[derive(Clone, Debug)]
struct GoCamNode {
    pub individual_gocam_id: IndividualId,
    pub id: String,
    pub label: String,
    pub node_type: GoCamNodeType,
    pub has_input: Vec<GoCamInput>,
    pub has_output: Vec<GoCamOutput>,
    pub located_in: Vec<GoCamComponent>,
    pub occurs_in: Vec<GoCamComponent>,
    pub part_of_process: Option<GoCamProcess>,
}

impl Display for GoCamNode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t", self.id)?;
        write!(f, "{}\t", self.label)?;
        write!(f, "{}\t", self.enabler_label())?;
        if let Some(ref part_of_process) = self.part_of_process {
            write!(f, "{}\t", part_of_process.label_or_id())?;
        } else {
            write!(f, "\t")?;
        }
        let has_input_string =
            self.has_input.iter().map(|l| l.label_or_id()).collect::<Vec<_>>().join(",");
        if has_input_string.len() > 0 {
            write!(f, "{}\t", has_input_string)?;
        } else {
            write!(f, "\t")?;
        }
        let has_output_string =
            self.has_output.iter().map(|l| l.label_or_id()).collect::<Vec<_>>().join(",");
        if has_output_string.len() > 0 {
            write!(f, "{}\t", has_output_string)?;
        } else {
            write!(f, "\t")?;
        }
        let occurs_in_string =
            self.occurs_in.iter().map(|l| l.label_or_id()).collect::<Vec<_>>().join(",");
        if occurs_in_string.len() > 0 {
            write!(f, "{}", occurs_in_string)?
        }
        let located_in_string =
            self.located_in.iter().map(|l| l.label_or_id()).collect::<Vec<_>>().join(",");
        if located_in_string.len() > 0 {
            write!(f, "{}\t", located_in_string)?;
        } else {
            write!(f, "\t")?;
        }
        Ok(())
    }
}

impl GoCamNode {
    pub fn enabler_label(&self) -> &str {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            enabler.label()
        } else {
            ""
        }
    }
}

#[derive(Clone, Debug)]
struct GoCamEdge {
    pub fact_id: FactId,
    pub relation_id: String,
    pub relation_label: String,
}

const MOLECULAR_FUNCTION_ID: &str = "GO:0003674";
const CELLULAR_COMPONENT_ID: &str = "GO:0032991";
const BIOLOGICAL_PROCESS_ID: &str = "GO:0008150";
const PROTEIN_CONTAINING_COMPLEX_ID: &str = "GO:0032991";
const CHEBI_PROTEIN_ID: &str = "CHEBI:36080";
const CHEBI_CHEMICAL_ENTITY_ID: &str = "CHEBI:24431";

fn has_root_term(individual: &Individual, term_id: &str) -> bool {
    for individual_type in &individual.root_types {
        if let Some(ref individual_type_id) = individual_type.id {
            if individual_type_id == term_id {
                return true;
            }
        }
    }

    false
}

fn individual_is_activity(individual: &Individual) -> bool {
    has_root_term(individual, MOLECULAR_FUNCTION_ID)
}

fn individual_is_component(individual: &Individual) -> bool {
    has_root_term(individual, CELLULAR_COMPONENT_ID)
}

fn individual_is_process(individual: &Individual) -> bool {
    has_root_term(individual, BIOLOGICAL_PROCESS_ID)
}

fn individual_is_complex(individual: &Individual) -> bool {
    has_root_term(individual, PROTEIN_CONTAINING_COMPLEX_ID)
}
*/

fn individual_is_chemical(individual: &Individual) -> bool {
    if !has_root_term(individual, CHEBI_CHEMICAL_ENTITY_ID) {
        return false;
    }

    if let Some(ref id) = get_individual_type(individual).id {
        if id.starts_with("CHEBI:") {
            return true;
        }
    }

    false
}

fn get_individual_type(individual: &Individual) -> &IndividualType {
    &individual.types[0]
}

fn individual_is_unknown_protein(individual: &Individual) -> bool {
    let individual_type = get_individual_type(individual);

    if let Some(ref individual_type_id) = individual_type.id {
        if individual_type_id == CHEBI_PROTEIN_ID {
            return true;
        }
    }

    false
}

fn is_gene_id(identifier: &str) -> bool {
    ["PomBase:", "FB:", "UniProtKB:"].iter().any(|s| identifier.starts_with(*s))
}

fn make_graph(model: &GoCamModel) -> GoCamGraph {
    let model_id = model.id();
    let model_title = model.title();

    let mut graph = GoCamGraph::new();

    let mut temp_nodes = HashMap::new();

    for individual in model.individuals() {
        if individual_is_activity(individual) ||
            individual_is_chemical(individual) &&
            !individual_is_unknown_protein(individual) {
                let individual_type = get_individual_type(individual);
                let detail =
                    if individual_is_chemical(individual) {
                        GoCamNodeType::Chemical
                    } else {
                        GoCamNodeType::Unknown
                    };
            let gocam_node = GoCamNode {
                individual_gocam_id: individual.id.clone(),
                id: individual_type.id.clone().unwrap_or_else(|| "NO_ID".to_owned()),
                label: individual_type.label.clone().unwrap_or_else(|| "NO_LABEL".to_owned()),
                node_type: detail,
                has_input: vec![],
                has_output: vec![],
                located_in: vec![],
                occurs_in: vec![],
                part_of_process: None,
            };

            temp_nodes.insert(individual.id.clone(), gocam_node);
        }
    }

    for fact in model.facts() {
        let Some(subject_node) = temp_nodes.get_mut(&fact.subject)
        else {
            continue;
        };

        let object_individual = model.fact_object(fact);
        let object_type = get_individual_type(object_individual);

        match fact.property_label.as_str() {
            "enabled by" => {
                if let Some(ref object_type_id) = object_type.id {
                    if is_gene_id(object_type_id) {
                        let gene_enabler = GoCamActivity::Gene(object_type.clone());
                        subject_node.node_type = GoCamNodeType::Activity(gene_enabler);
                    }
                    else if object_type_id.starts_with("CHEBI:") {
                        let chemical_enabler = GoCamActivity::Chemical(object_type.clone());
                        subject_node.node_type = GoCamNodeType::Activity(chemical_enabler);
                    }
                    else if object_type_id.starts_with("GO:") || object_type_id.starts_with("ComplexPortal:") {
                        let complex_enabler = GoCamActivity::Complex(object_type.clone());
                        subject_node.node_type = GoCamNodeType::Activity(complex_enabler);
                    }
                    else if object_type_id.starts_with("PR:") {
                        let modified_protein_enabler = GoCamActivity::ModifiedProtein(object_type.clone());
                        subject_node.node_type = GoCamNodeType::Activity(modified_protein_enabler);
                    }
                    else  {
                        eprintln!("can't handle enabled by object: {}", object_individual.id);
                    }
                }
            },
            "has input" => {
                subject_node.has_input.push(object_type.clone());
            },
            "has output" => {
                subject_node.has_output.push(object_type.clone());
            },
            "located in" => {
                subject_node.located_in.push(object_type.clone());
            },
            "occurs in" => {
                subject_node.occurs_in.push(object_type.clone());
            },
            "part of" => {
                subject_node.part_of_process = Some(object_type.clone());
            },
            &_ => {
                // eprintln!("ignoring rel from fact: {} {}", fact.property_label, fact.id());
            }
        }
    }

    let mut id_map = HashMap::new();

    for node in temp_nodes.values() {
        let idx = graph.add_node(node.to_owned());
        id_map.insert(node.individual_gocam_id.clone(), idx);
    }

    for fact in model.facts() {
        let subject_id = &fact.subject;
        let object_id = &fact.object;

        let subject_node = temp_nodes.get(subject_id);
        let object_node = temp_nodes.get(object_id);

        if let (Some(subject_node), Some(object_node)) =
              (temp_nodes.get(subject_id), temp_nodes.get(object_id))
        {
/*
            println!("{} ({}) <- {} -> {} ({})",
                     subject_node.label(), subject_node.individual_type.label_or_id(),
                     fact.property_label,
                object_node.label(), object_node.individual_type.label_or_id());
 
*/
        }
    }

    for node in temp_nodes.values() {
        if node.enabler_label() == "protein" {
            println!("{}\t{}\t{}\t{}", model_id, model_title, node.individual_gocam_id, node);
        }
    }

    graph
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

fn print_stats(_model: &GoCamModel) {

}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    match args.action {
        Action::Stats { paths } => {
            for path in paths {
                let mut source = File::open(path).unwrap();
                let model = gocam_parse(&mut source)?;
                println!("{}", model.id());

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
        }
    }

    Ok(())
}
