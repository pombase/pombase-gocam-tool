use pombase_gocam::GoCamModel;

use crate::{allowed_relation_config::AllowedRelConfig, ontology_closure::OntologyClosure};

use pombase_gocam::REL_NAMES;

pub(crate) fn check_relations(model: &GoCamModel,
                              config: &AllowedRelConfig,
                              ontology_closure: &OntologyClosure)
    -> Vec<String>
{
    let mut ret = vec![];

    for (_, source_idx, edge, object_idx) in model.edge_iterator() {
        let subject_node = model.graph().node_weight(source_idx).unwrap();
        let object_node = model.graph().node_weight(object_idx).unwrap();
        if !subject_node.is_activity() || !object_node.is_activity() {
            continue;
        }
        let subject_term_id = &subject_node.node_id;
        let rel_id = &edge.id;
        let Some(term_id_and_parents) = ontology_closure.get_term_parents(subject_term_id)
        else {
            ret.push(format!("can't find parents of {}", subject_term_id));
            continue;
        };

        for (parent_rel, parent_term_id) in term_id_and_parents {
            if parent_rel != "rdfs:subClassOf" {
                continue;
            }

            let Some(allowed_term_rels) = config.term_config().get(parent_term_id.as_str())
            else {
                ret.push(format!("no rels for {} ({})", parent_term_id, subject_term_id));
                continue;
            };

            if !allowed_term_rels.contains(rel_id) {
                let rel_name = REL_NAMES.get(rel_id).unwrap();
                ret.push(format!("relation {} ({}) not allowed for {} because of config for {}",
                         rel_name, rel_id, subject_term_id, parent_term_id));
            }
        }
    }

    ret
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use std::fs::File;

    use pombase_gocam::parse_gocam_py_model;

    use super::check_relations;

    #[test]
    fn check_test() {
        let closure_file = File::open("tests/data/closure.tsv").unwrap();
        let mut closure_buf_reader = BufReader::new(closure_file);
        let closure = crate::parse_closure(&mut closure_buf_reader).unwrap();

        let config_file = File::open("tests/data/config.tsv").unwrap();
        let mut config_buf_reader = BufReader::new(config_file);
        let config = crate::parse_allowed_relations_config(&mut config_buf_reader).unwrap();

        let mut source = File::open("tests/data/67ae98b500000055.yaml").unwrap();
        let model = parse_gocam_py_model(&mut source).unwrap();

        let warnings = check_relations(&model, &config, &closure);

        assert_eq!(warnings.len(), 2);
        assert_eq!(warnings[0], "relation directly negatively regulates (RO:0002630) not allowed for GO:0003713 because of config for GO:0003824");
    }
}
