use std::{collections::{HashMap, HashSet}, io::BufRead};

use thiserror::Error;

#[derive(Error, Debug)]
pub enum ClosureError {
    #[error("line doesn't have 3 fields")]
    FieldCountError { detail: String },
    #[error("I/O error: {0}")]
    IOError(#[from] std::io::Error),
}

pub type TermId = String;
pub type RelId = String;

pub(crate) struct OntologyClosure {
    pub term_parents: HashMap<TermId, HashSet<(RelId, TermId)>>,
}

pub(crate) fn parse_closure(buf_reader: &mut dyn BufRead)
    -> Result<OntologyClosure, ClosureError>
{
    let mut term_parents = HashMap::new();

    let mut all_term_ids = HashSet::new();

    for line_result in buf_reader.lines() {
        let line = line_result?;

        let bits: Vec<_> = line.split("\t").collect();

        if bits.len() != 3 {
            let detail = format!("line doesn't have 3 fields: {}", line);
            return Err(ClosureError::FieldCountError {
                detail
            });
        }

        let subject = bits[0];
        all_term_ids.insert(subject.to_owned());
        let rel = bits[1].to_owned();
        let object = bits[2];
        all_term_ids.insert(object.to_owned());

        term_parents.entry(subject.to_owned())
            .or_insert_with(HashSet::new)
            .insert((rel, object.to_owned()));
    }

    for term_id in all_term_ids.into_iter() {
        term_parents.entry(term_id.clone())
            .or_insert_with(HashSet::new)
            .insert(("rdfs:subClassOf".to_owned(), term_id));
    }

    Ok(OntologyClosure {
        term_parents,
    })
}

impl OntologyClosure {
    pub fn get_term_parents(&self, id: &str) -> Option<&HashSet<(RelId, TermId)>> {
        self.term_parents.get(id)
    }
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use std::fs::File;

    #[test]
    fn parse_test() {
        let file = File::open("tests/data/closure.tsv").unwrap();
        let mut buf_reader = BufReader::new(file);
        let closure = crate::parse_closure(&mut buf_reader).unwrap();

        assert_eq!(closure.term_parents.len(), 18);

        let go_0070858_parents = closure.term_parents.get("GO:0070858").unwrap();
        assert_eq!(go_0070858_parents.len(), 1);

        assert_eq!(go_0070858_parents.iter().next().unwrap(),
                   &("RO:0002212".to_owned(), "GO:0032787".to_owned()));

        let go_1901953_parents = closure.get_term_parents("GO:1901953").unwrap();
        assert_eq!(go_1901953_parents.len(), 2);
    }
}
