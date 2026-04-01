use std::{collections::HashMap, io::{BufWriter, Error, Write}};

use chrono::{DateTime, Local};
use itertools::Itertools;
use pombase_gocam::gocam_py::{EvidenceItem, GoCamPyModel};

fn make_annotation_line(db_name: &str, gene_uniquename: &str,
                        term_id: &str, aspect: &str,
                        reference: &str, ev_code: &str, with_from: &str,
                        taxon_id: &str, date: &str)
    -> String
{
    format!("{}\t{}\t\t\t{}\t{}\t{}\t{}\t{}\t\t\tprotein\t{}\t{}\t{}\t\t",
            db_name, gene_uniquename, term_id, reference, ev_code, with_from, aspect, taxon_id, date, db_name)
}


pub fn write_go_annotation_file(writer: &mut dyn Write,
                                go_evidence_code_map: &HashMap<String, String>,
                                model: &GoCamPyModel, db_name: &str)
   -> Result<(), Error>
{
    let mut buf_writer = BufWriter::new(writer);

    let Some(ref taxon_id) = model.taxon
    else {
        return Ok(());
    };

    let taxon_id = taxon_id.replace("NCBITaxon:", "taxon:");

    let date_modified = if let Some(ref date_modified) = model.date_modified {
        date_modified.to_owned()
    } else {
        let local: DateTime<Local> = Local::now();
        local.format("%F").to_string()
    };

    let details_from_item = |evidence_item: &EvidenceItem| {
        let Some(reference) = &evidence_item.reference
        else {
            return None;
        };

        let ev_code = &evidence_item.term;

        if !reference.starts_with("GO_REF:") {
            return None;
        }

        let Some(go_ev_code) = go_evidence_code_map.get(ev_code)
        else {
            panic!("unknown evidence code: {}", ev_code);
        };

        let with_from = evidence_item.with_objects.iter().join(",");

        Some((go_ev_code, reference.to_owned(), with_from))
    };


    for activity in &model.activities {
        let gene_uniquename = &activity.enabled_by.term;
        let mf_term_id = &activity.molecular_function.term;

        if let Some(evidence_item) = activity.enabled_by.evidence.get(0) {
            let Some((go_ev_code, reference, with_from)) = details_from_item(evidence_item)
            else {
                continue;
            };

            let line = make_annotation_line(db_name, gene_uniquename, mf_term_id,
                                            "F", &reference,
                                            go_ev_code, &with_from, &taxon_id, &date_modified);

            writeln!(buf_writer, "{}", line)?;
        }

        if let Some(ref bp_association) = activity.part_of &&
           let Some(evidence_item) = bp_association.evidence.get(0) {

            let bp_term_id = &bp_association.term;

            let Some((go_ev_code, reference, with_from)) = details_from_item(evidence_item)
            else {
                continue;
            };

            let line = make_annotation_line(db_name, gene_uniquename, bp_term_id,
                                            "P", &reference,
                                            go_ev_code, &with_from, &taxon_id, &date_modified);

            writeln!(buf_writer, "{}", line)?;

        }

        if let Some(ref cc_assocation) = activity.occurs_in &&
            let Some(evidence_item) = cc_assocation.evidence.get(0) {

            let cc_term_id = &cc_assocation.term;

            let Some((go_ev_code, reference, with_from)) = details_from_item(evidence_item)
            else {
                continue;
            };

            let line = make_annotation_line(db_name, gene_uniquename, cc_term_id,
                                            "C", &reference,
                                            go_ev_code, &with_from, &taxon_id, &date_modified);

            writeln!(buf_writer, "{}", line)?;
        }
    }

    Ok(())
}
