use std::collections::HashSet;

use pombase_gocam::GoCamModel;

use crate::ontology_info::OntologyInfo;

pub type TermId = String;
pub type TermName = String;
pub type FieldType = String;

pub(crate) fn find_obsolete_terms(model: &GoCamModel,
                                  ontology_info: &OntologyInfo)
    -> Vec<(TermId, TermName, FieldType, Option<String>)>
{
    let mut seen = HashSet::new();
    let mut ret = vec![];

    for (_, node) in model.node_iterator() {
        if node.is_activity() {
            if ontology_info.is_obsolete_term(&node.node_id) {
                let obsolete = (node.node_id.clone(), node.label.clone(), "molecular_function".to_owned(), None);
                if !seen.contains(&obsolete) {
                    ret.push(obsolete.clone());
                    seen.insert(obsolete);
                }
            }

            for occurs_in in node.occurs_in.iter() {
                if ontology_info.is_obsolete_term(occurs_in.id()) {
                    let obsolete = (occurs_in.id().to_owned(), occurs_in.label().to_owned(), "occurs_in".to_owned(), None);
                    if !seen.contains(&obsolete) {
                        ret.push(obsolete.clone());
                        seen.insert(obsolete);
                    }
                }
            }

            if let Some(ref part_of_process) = node.part_of_process &&
                ontology_info.is_obsolete_term(part_of_process.id())
            {
                let obsolete = (part_of_process.id().to_owned(), part_of_process.label.to_owned(),
                                "part_of_process".to_owned(), None);
                if !seen.contains(&obsolete) {
                    ret.push(obsolete.clone());
                    seen.insert(obsolete);
                }
            }

        }
    }

    ret
}
