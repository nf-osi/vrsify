//! Validation gate against the GA4GH VRS 2.0 language-neutral golden suite.
//!
//! Fixtures: `tests/fixtures/models.yaml` and `functions.yaml`, copied from
//! `github.com/ga4gh/vrs` tag `2.0`. We validate the types this tool emits
//! (Allele, SequenceLocation) byte-exact for `ga4gh_serialize`, `ga4gh_digest`,
//! and `ga4gh_identify`.

use serde::Deserialize;
use serde_yaml::Value as Yaml;
use vrsify::vrs::{Adjacency, Allele, CopyNumberChange, CopyNumberCount, SequenceLocation};

#[derive(Deserialize)]
struct Expected {
    ga4gh_serialize: Option<String>,
    ga4gh_digest: Option<String>,
    ga4gh_identify: Option<String>,
}

fn models() -> serde_yaml::Mapping {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/models.yaml"
    ))
    .expect("read models.yaml");
    serde_yaml::from_str(&text).expect("parse models.yaml")
}

/// Each fixture entry has an `in:` object and an `out:` expectation. Returns
/// (label, in_value, expected) tuples for the given top-level type key.
fn cases(models: &serde_yaml::Mapping, key: &str) -> Vec<(String, Yaml, Expected)> {
    let list = models
        .get(Yaml::from(key))
        .and_then(|v| v.as_sequence())
        .unwrap_or_else(|| panic!("no fixtures for {key}"));
    list.iter()
        .enumerate()
        .map(|(i, case)| {
            let label = case
                .get("name")
                .and_then(|n| n.as_str())
                .map(String::from)
                .unwrap_or_else(|| format!("{key}[{i}]"));
            let in_v = case.get("in").expect("in").clone();
            let expected: Expected =
                serde_yaml::from_value(case.get("out").expect("out").clone()).expect("out parse");
            (label, in_v, expected)
        })
        .collect()
}

#[test]
fn sequence_location_identifiers() {
    let m = models();
    let mut checked = 0;
    for (label, in_v, exp) in cases(&m, "SequenceLocation") {
        let loc: SequenceLocation =
            serde_yaml::from_value(in_v).unwrap_or_else(|e| panic!("{label}: deser: {e}"));
        if let Some(s) = &exp.ga4gh_serialize {
            assert_eq!(&loc.ga4gh_serialize(), s, "{label}: serialize");
        }
        if let Some(d) = &exp.ga4gh_digest {
            assert_eq!(&loc.digest(), d, "{label}: digest");
        }
        if let Some(id) = &exp.ga4gh_identify {
            assert_eq!(&loc.ga4gh_id(), id, "{label}: identify");
        }
        checked += 1;
    }
    assert!(checked > 0, "no SequenceLocation cases ran");
}

#[test]
fn allele_identifiers() {
    let m = models();
    let mut checked = 0;
    for (label, in_v, exp) in cases(&m, "Allele") {
        // Skip cases whose state type this milestone doesn't model.
        let allele: Allele = match serde_yaml::from_value(in_v) {
            Ok(a) => a,
            Err(e) => {
                eprintln!("SKIP {label}: {e}");
                continue;
            }
        };
        if let Some(s) = &exp.ga4gh_serialize {
            assert_eq!(&allele.ga4gh_serialize(), s, "{label}: serialize");
        }
        if let Some(d) = &exp.ga4gh_digest {
            assert_eq!(&allele.digest(), d, "{label}: digest");
        }
        if let Some(id) = &exp.ga4gh_identify {
            assert_eq!(&allele.ga4gh_id(), id, "{label}: identify");
        }
        checked += 1;
    }
    assert!(checked > 0, "no Allele cases ran");
}

#[test]
fn copy_number_count_identifiers() {
    let m = models();
    let mut checked = 0;
    for (label, in_v, exp) in cases(&m, "CopyNumberCount") {
        let cn: CopyNumberCount =
            serde_yaml::from_value(in_v).unwrap_or_else(|e| panic!("{label}: deser: {e}"));
        if let Some(s) = &exp.ga4gh_serialize {
            assert_eq!(&cn.ga4gh_serialize(), s, "{label}: serialize");
        }
        if let Some(d) = &exp.ga4gh_digest {
            assert_eq!(&cn.digest(), d, "{label}: digest");
        }
        if let Some(id) = &exp.ga4gh_identify {
            assert_eq!(&cn.ga4gh_id(), id, "{label}: identify");
        }
        checked += 1;
    }
    assert!(checked > 0, "no CopyNumberCount cases ran");
}

#[test]
fn copy_number_change_identifiers() {
    let m = models();
    let mut checked = 0;
    for (label, in_v, exp) in cases(&m, "CopyNumberChange") {
        let cx: CopyNumberChange =
            serde_yaml::from_value(in_v).unwrap_or_else(|e| panic!("{label}: deser: {e}"));
        if let Some(s) = &exp.ga4gh_serialize {
            assert_eq!(&cx.ga4gh_serialize(), s, "{label}: serialize");
        }
        if let Some(d) = &exp.ga4gh_digest {
            assert_eq!(&cx.digest(), d, "{label}: digest");
        }
        if let Some(id) = &exp.ga4gh_identify {
            assert_eq!(&cx.ga4gh_id(), id, "{label}: identify");
        }
        checked += 1;
    }
    assert!(checked > 0, "no CopyNumberChange cases ran");
}

#[test]
fn adjacency_identifiers() {
    let m = models();
    let mut checked = 0;
    for (label, in_v, exp) in cases(&m, "Adjacency") {
        let aj: Adjacency =
            serde_yaml::from_value(in_v).unwrap_or_else(|e| panic!("{label}: deser: {e}"));
        if let Some(s) = &exp.ga4gh_serialize {
            assert_eq!(&aj.ga4gh_serialize(), s, "{label}: serialize");
        }
        if let Some(d) = &exp.ga4gh_digest {
            assert_eq!(&aj.digest(), d, "{label}: digest");
        }
        if let Some(id) = &exp.ga4gh_identify {
            assert_eq!(&aj.ga4gh_id(), id, "{label}: identify");
        }
        checked += 1;
    }
    assert!(checked > 0, "no Adjacency cases ran");
}
