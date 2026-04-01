use std::{collections::HashMap, io::{BufWriter, Error, Write}};

use chrono::{DateTime, Local};
use itertools::Itertools;
use pombase_gocam::gocam_py::GoCamPyModel;

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

    for activity in &model.activities {
        let gene_uniquename = &activity.enabled_by.term;
        let mf_term_id = &activity.molecular_function.term;

        if let Some(mf_evidence_item) = activity.enabled_by.evidence.get(0) &&
           let Some(reference) = &mf_evidence_item.reference {

            let ev_code = &mf_evidence_item.term;
            let with_from = mf_evidence_item.with_objects.iter().join(",");

            if !reference.starts_with("GO_REF:") {
                continue;
            }

            let Some(go_ev_code) = go_evidence_code_map.get(ev_code)
            else {
                panic!("unknown evidence code: {}", ev_code);
            };

            let line = make_annotation_line(db_name, gene_uniquename, mf_term_id,
                                            "F", reference,
                                            go_ev_code, &with_from, &taxon_id, &date_modified);

            writeln!(buf_writer, "{}", line)?;
        }
    }

    Ok(())
}
