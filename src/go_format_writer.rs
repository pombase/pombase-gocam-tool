use std::{collections::HashMap, io::{BufWriter, Error, Write}};

use chrono::{DateTime, Local};
use itertools::Itertools;
use pombase_gocam::gocam_py::{BiologicalProcessAssociation, EvidenceItem, GoCamPyModel};

const ND_ECO_TERM_ID: &str = "ECO:0000307";

#[allow(clippy::too_many_arguments)]
fn make_annotation_line(db_name: &str, gene_uniquename: &str,
                        term_id: &str, aspect: &str,
                        reference: &str, ev_code: &str, with_from: &str,
                        taxon_id: &str, date: &str)
    -> String
{
    format!("{}\t{}\t\t\t{}\t{}\t{}\t{}\t{}\t\t\tprotein\t{}\t{}\t{}\t\t",
            db_name, gene_uniquename, term_id, reference, ev_code, with_from, aspect, taxon_id, date, db_name)
}

fn details_from_item(go_evidence_code_map: &HashMap<String, String>,
                     evidence_item: &EvidenceItem)
    -> Result<Option<(String, String, String)>, String>
{
    let Some(reference) = &evidence_item.reference
    else {
        return Ok(None);
    };

    let ev_code = &evidence_item.term;

    if ev_code == ND_ECO_TERM_ID {
        return Ok(None);
    }

    let Some(go_ev_code) = go_evidence_code_map.get(ev_code)
    else {
        return Err(format!("unknown evidence code: {}", ev_code));
    };

    if !reference.starts_with("GO_REF:") && ev_code != "ECO:0000304" {
        return Ok(None);
    }

    let with_from = evidence_item.with_objects.iter().join(",");

    Ok(Some((go_ev_code.to_owned(), reference.to_owned(), with_from)))
}


fn bp_association_lines(go_ev_code_map: &HashMap<String, String>,
                        db_name: &str, taxon_id: &str,
                        gene_uniquename: &str, date_modified: &str,
                        bp_association: &BiologicalProcessAssociation)
    -> Result<Vec<String>, String>
{
    let mut ret = vec![];

    if let Some(evidence_item) = bp_association.evidence.first()
    {

        let bp_term_id = &bp_association.term;

        let details = details_from_item(go_ev_code_map, evidence_item)?;

        if let Some((go_ev_code, reference, with_from)) = details {
            let line = make_annotation_line(db_name, gene_uniquename, bp_term_id,
                                            "P", &reference,
                                            &go_ev_code, &with_from, taxon_id, date_modified);

            ret.push(line)
        }

        if let Some(ref bp_association) = bp_association.part_of {
            let lines =
                bp_association_lines(go_ev_code_map, db_name, taxon_id,
                                     gene_uniquename, date_modified, bp_association)?;
            ret.extend_from_slice(&lines);
        }
    }

    Ok(ret)
}

pub fn write_go_annotation_file(writer: &mut dyn Write,
                                go_ev_code_map: &HashMap<String, String>,
                                model: &GoCamPyModel, db_name: &str)
   -> Result<(), Error>
{
    let mut buf_writer = BufWriter::new(writer);

    let Some(ref taxon_id) = model.taxon
    else {
        return Ok(());
    };

    let taxon_id = taxon_id.replace("NCBITaxon:", "taxon:");

    let db_name_prefix = db_name.to_owned() + ":";

    let date_modified = if let Some(ref date_modified) = model.date_modified {
        date_modified.to_owned()
    } else {
        let local: DateTime<Local> = Local::now();
        local.format("%F").to_string()
    };

    for activity in &model.activities {
        let Some(gene_uniquename) = &activity.enabled_by.term.strip_prefix(&db_name_prefix)
        else {
            // not a gene
            continue;
        };

        let mf_term_id = &activity.molecular_function.term;

        if let Some(evidence_item) = activity.enabled_by.evidence.first() {
            match details_from_item(go_ev_code_map, evidence_item) {
                Ok(details) => {
                    if let Some((go_ev_code, reference, with_from)) = details {
                        let line = make_annotation_line(db_name, gene_uniquename, mf_term_id,
                                                        "F", &reference,
                                                        &go_ev_code, &with_from, &taxon_id, &date_modified);

                        writeln!(buf_writer, "{}", line)?;
                    }
                },
                Err(message) => {
                    eprintln!("warning in {} for {} {}: {}",
                              model.id, mf_term_id, gene_uniquename, message);
                }
            };
        }

        if let Some(ref bp_association) = activity.part_of {
            let bp_lines_result =
                bp_association_lines(go_ev_code_map, db_name, &taxon_id,
                                     gene_uniquename, &date_modified, bp_association);

            match bp_lines_result {
                Ok(bp_lines) => {
                    for line in bp_lines {
                        writeln!(buf_writer, "{}", line)?;
                    }
                },
                Err(message) => {
                    eprintln!("warning in {} for {} {}: {}",
                              model.id, mf_term_id, gene_uniquename, message);
                }
            };
        }

        if let Some(ref cc_assocation) = activity.occurs_in &&
            let Some(evidence_item) = cc_assocation.evidence.first() {

            let cc_term_id = &cc_assocation.term;

            let details = match details_from_item(go_ev_code_map, evidence_item) {
                Ok(details) => details,
                Err(message) => {
                    eprintln!("warning in {} for {} {}: {}",
                              model.id, mf_term_id, gene_uniquename, message);
                    continue;
                }
            };

            let Some((go_ev_code, reference, with_from)) = details
            else {
                continue;
            };

            let line = make_annotation_line(db_name, gene_uniquename, cc_term_id,
                                            "C", &reference,
                                            &go_ev_code, &with_from, &taxon_id, &date_modified);

            writeln!(buf_writer, "{}", line)?;
        }
    }

    Ok(())
}
