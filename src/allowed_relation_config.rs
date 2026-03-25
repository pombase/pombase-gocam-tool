use std::{collections::{HashMap, HashSet}, io::BufRead};

use thiserror::Error;

#[derive(Error, Debug)]
pub enum AllowedRelConfigError {
    #[error("unknown relation name")]
    UnknownRelationNameError { detail: String },
    #[error("I/O error: {0}")]
    IOError(#[from] std::io::Error),
}

use pombase_gocam::REL_NAMES;

pub type TermId = String;
pub type RelId = String;

pub(crate) struct AllowedRelConfig {
    pub _term_config: HashMap<TermId, HashSet<RelId>>,
}

impl AllowedRelConfig {
    pub fn term_config(&self) -> &HashMap<TermId, HashSet<RelId>> {
        &self._term_config
    }
}

pub(crate) fn parse_allowed_relations_config(buf_reader: &mut dyn BufRead)
    -> Result<AllowedRelConfig, AllowedRelConfigError>
{
    let rel_ids_by_name: HashMap<String, RelId> =
        REL_NAMES.entries()
        .map(|(k, v)| (v.to_string(), k.to_string())).collect();

    let mut term_config = HashMap::new();

    for line_result in buf_reader.lines() {
        let line = line_result?;

        let bits: Vec<_> = line.split("\t").collect();

        let term_id = bits[0];

        if term_id == "term" {
            // header
            continue;
        }

        let allowed_rel_names = bits[2].to_owned();

        for allowed_rel_name in allowed_rel_names.split("|") {
            let Some(rel_id) = rel_ids_by_name.get(&allowed_rel_name.replace("_", " "))
            else {
                let detail = format!("unknown relation name: {}", allowed_rel_name);
                return Err(AllowedRelConfigError::UnknownRelationNameError {
                    detail
                });

            };
            term_config.entry(term_id.to_owned())
                .or_insert_with(HashSet::new)
                .insert(rel_id.to_owned());
        }
    }

    Ok(AllowedRelConfig {
        _term_config: term_config,
    })
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use std::io::BufReader;
    use std::fs::File;

    #[test]
    fn parse_test() {
        let file = File::open("tests/data/config.tsv").unwrap();
        let mut buf_reader = BufReader::new(file);
        let config = crate::parse_allowed_relations_config(&mut buf_reader).unwrap();

        assert_eq!(config.term_config().len(), 3);

        let rels = config.term_config().get("GO:0003824").unwrap();

        let mut expected_rels = HashSet::new();
        expected_rels.insert("RO:0012009".to_owned());
        expected_rels.insert("RO:0002629".to_owned());
        expected_rels.insert("RO:0002409".to_owned());
        expected_rels.insert("RO:0002413".to_owned());
        assert_eq!(rels, &expected_rels);
    }
}
