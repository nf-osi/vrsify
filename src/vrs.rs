//! Minimal GA4GH VRS 2.0 model with computed-identifier support.
//!
//! Covers the types needed for SNV/indel ingestion from VCF: `Allele`,
//! `SequenceLocation`, `SequenceReference`, and the `LiteralSequenceExpression` /
//! `ReferenceLengthExpression` / `LengthExpression` states. It also models the
//! copy-number branch (`CopyNumberCount` / `CopyNumberChange`) and the structural
//! `Adjacency` type, which back CNV/SV ingestion.
//!
//! The `ga4gh_serialize` / digest / identify logic is validated byte-exact against the
//! GA4GH `vrs/validation/models.yaml` golden fixtures (see `tests/`).

use serde::Deserialize;
use serde_json::{Map, Value};

use crate::digest::sha512t24u;

/// A VRS coordinate: either a definite integer or a `[min, max]` range whose ends may
/// be indefinite (`null`). Interbase (0-based) throughout.
#[derive(Debug, Clone, Deserialize)]
#[serde(untagged)]
pub enum Coordinate {
    Definite(i64),
    Range([Option<i64>; 2]),
}

impl Coordinate {
    fn to_value(&self) -> Value {
        match self {
            Coordinate::Definite(n) => Value::from(*n),
            Coordinate::Range([a, b]) => Value::Array(vec![opt(a), opt(b)]),
        }
    }
}

fn opt(v: &Option<i64>) -> Value {
    match v {
        Some(n) => Value::from(*n),
        None => Value::Null,
    }
}

/// A reference to a sequence by its refget accession (`SQ.<digest>`). Non-identifiable:
/// it has no computed digest of its own and is inlined wherever it appears.
#[derive(Debug, Clone, Deserialize)]
pub struct SequenceReference {
    #[serde(rename = "refgetAccession")]
    pub refget_accession: String,
}

impl SequenceReference {
    fn digest_value(&self) -> Value {
        let mut m = Map::new();
        m.insert("type".into(), Value::from("SequenceReference"));
        m.insert(
            "refgetAccession".into(),
            Value::from(self.refget_accession.clone()),
        );
        Value::Object(m)
    }
}

/// The allelic state. Non-identifiable; inlined into the Allele.
#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "type")]
pub enum State {
    LiteralSequenceExpression {
        sequence: String,
    },
    ReferenceLengthExpression {
        length: i64,
        #[serde(rename = "repeatSubunitLength")]
        repeat_subunit_length: i64,
        /// The reference-derived alt sequence. Present for output fidelity but
        /// deliberately EXCLUDED from the digest (see `digest_value`), matching VRS 2.0.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        sequence: Option<String>,
    },
    LengthExpression {
        length: i64,
    },
}

impl State {
    fn digest_value(&self) -> Value {
        let mut m = Map::new();
        match self {
            State::LiteralSequenceExpression { sequence } => {
                m.insert("type".into(), Value::from("LiteralSequenceExpression"));
                m.insert("sequence".into(), Value::from(sequence.clone()));
            }
            // Note: RLE's `sequence` is deliberately excluded from the digest.
            State::ReferenceLengthExpression {
                length,
                repeat_subunit_length,
                ..
            } => {
                m.insert("type".into(), Value::from("ReferenceLengthExpression"));
                m.insert("length".into(), Value::from(*length));
                m.insert("repeatSubunitLength".into(), Value::from(*repeat_subunit_length));
            }
            State::LengthExpression { length } => {
                m.insert("type".into(), Value::from("LengthExpression"));
                m.insert("length".into(), Value::from(*length));
            }
        }
        Value::Object(m)
    }

    /// Full JSON for NDJSON output. Identical to `digest_value` except that an RLE's
    /// reference-derived `sequence` is re-attached (it is present in output but excluded
    /// from the digest per VRS 2.0), matching vrs-python's default `rle_seq_limit=50`.
    fn output_value(&self) -> Value {
        let mut v = self.digest_value();
        if let State::ReferenceLengthExpression {
            sequence: Some(s), ..
        } = self
            && let Value::Object(m) = &mut v {
                m.insert("sequence".into(), Value::from(s.clone()));
            }
        v
    }
}

/// A location on a reference sequence. Identifiable — prefix `SL`.
#[derive(Debug, Clone, Deserialize)]
pub struct SequenceLocation {
    #[serde(rename = "sequenceReference")]
    pub sequence_reference: SequenceReference,
    #[serde(default)]
    pub start: Option<Coordinate>,
    #[serde(default)]
    pub end: Option<Coordinate>,
}

impl SequenceLocation {
    fn digest_value(&self) -> Value {
        // `start` and `end` are *inherent* keys for SequenceLocation and are always
        // emitted into the digest serialization — `null` when a coordinate is absent
        // (e.g. one-sided locations used by Adjacency). Objects that supply both (the
        // SNV/indel path) are unaffected.
        let mut m = Map::new();
        m.insert(
            "end".into(),
            self.end.as_ref().map(|e| e.to_value()).unwrap_or(Value::Null),
        );
        m.insert(
            "sequenceReference".into(),
            self.sequence_reference.digest_value(),
        );
        m.insert(
            "start".into(),
            self.start.as_ref().map(|s| s.to_value()).unwrap_or(Value::Null),
        );
        m.insert("type".into(), Value::from("SequenceLocation"));
        Value::Object(m)
    }

    pub fn ga4gh_serialize(&self) -> String {
        serde_json::to_string(&self.digest_value()).expect("serialize SequenceLocation")
    }

    pub fn digest(&self) -> String {
        sha512t24u(self.ga4gh_serialize().as_bytes())
    }

    pub fn ga4gh_id(&self) -> String {
        format!("ga4gh:SL.{}", self.digest())
    }

    /// Full VRS 2.0 JSON object with computed identifier, for embedding in output.
    pub fn to_output_value(&self) -> Value {
        let mut m = Map::new();
        m.insert("id".into(), Value::from(self.ga4gh_id()));
        m.insert("type".into(), Value::from("SequenceLocation"));
        if let Some(s) = &self.start {
            m.insert("start".into(), s.to_value());
        }
        if let Some(e) = &self.end {
            m.insert("end".into(), e.to_value());
        }
        m.insert(
            "sequenceReference".into(),
            self.sequence_reference.digest_value(),
        );
        Value::Object(m)
    }
}

/// A contextual, single molecular variation. Identifiable — prefix `VA`.
#[derive(Debug, Clone, Deserialize)]
pub struct Allele {
    pub location: SequenceLocation,
    pub state: State,
}

impl Allele {
    fn digest_value(&self) -> Value {
        let mut m = Map::new();
        m.insert("type".into(), Value::from("Allele"));
        // A nested identifiable object is referenced by its bare digest.
        m.insert("location".into(), Value::from(self.location.digest()));
        m.insert("state".into(), self.state.digest_value());
        Value::Object(m)
    }

    pub fn ga4gh_serialize(&self) -> String {
        serde_json::to_string(&self.digest_value()).expect("serialize Allele")
    }

    pub fn digest(&self) -> String {
        sha512t24u(self.ga4gh_serialize().as_bytes())
    }

    pub fn ga4gh_id(&self) -> String {
        format!("ga4gh:VA.{}", self.digest())
    }

    /// Full VRS 2.0 JSON object with computed identifiers, for NDJSON output.
    pub fn to_output_value(&self) -> Value {
        let mut m = Map::new();
        m.insert("id".into(), Value::from(self.ga4gh_id()));
        m.insert("type".into(), Value::from("Allele"));
        m.insert("location".into(), self.location.to_output_value());
        m.insert("state".into(), self.state.output_value());
        Value::Object(m)
    }
}

/// A copy-number count variation — an absolute count (or range) of copies of a
/// region. Identifiable — prefix `CN`.
#[derive(Debug, Clone, Deserialize)]
pub struct CopyNumberCount {
    pub location: SequenceLocation,
    pub copies: Coordinate,
}

impl CopyNumberCount {
    fn digest_value(&self) -> Value {
        let mut m = Map::new();
        m.insert("copies".into(), self.copies.to_value());
        m.insert("location".into(), Value::from(self.location.digest()));
        m.insert("type".into(), Value::from("CopyNumberCount"));
        Value::Object(m)
    }

    pub fn ga4gh_serialize(&self) -> String {
        serde_json::to_string(&self.digest_value()).expect("serialize CopyNumberCount")
    }

    pub fn digest(&self) -> String {
        sha512t24u(self.ga4gh_serialize().as_bytes())
    }

    pub fn ga4gh_id(&self) -> String {
        format!("ga4gh:CN.{}", self.digest())
    }

    /// Full VRS 2.0 JSON object with computed identifiers, for NDJSON output.
    pub fn to_output_value(&self) -> Value {
        let mut m = Map::new();
        m.insert("id".into(), Value::from(self.ga4gh_id()));
        m.insert("type".into(), Value::from("CopyNumberCount"));
        m.insert("location".into(), self.location.to_output_value());
        m.insert("copies".into(), self.copies.to_value());
        Value::Object(m)
    }
}

/// A copy-number change variation — a qualitative change relative to reference.
/// Identifiable — prefix `CX`.
#[derive(Debug, Clone, Deserialize)]
pub struct CopyNumberChange {
    pub location: SequenceLocation,
    #[serde(rename = "copyChange")]
    pub copy_change: String,
}

impl CopyNumberChange {
    fn digest_value(&self) -> Value {
        let mut m = Map::new();
        m.insert("copyChange".into(), Value::from(self.copy_change.clone()));
        m.insert("location".into(), Value::from(self.location.digest()));
        m.insert("type".into(), Value::from("CopyNumberChange"));
        Value::Object(m)
    }

    pub fn ga4gh_serialize(&self) -> String {
        serde_json::to_string(&self.digest_value()).expect("serialize CopyNumberChange")
    }

    pub fn digest(&self) -> String {
        sha512t24u(self.ga4gh_serialize().as_bytes())
    }

    pub fn ga4gh_id(&self) -> String {
        format!("ga4gh:CX.{}", self.digest())
    }

    /// Full VRS 2.0 JSON object with computed identifiers, for NDJSON output.
    pub fn to_output_value(&self) -> Value {
        let mut m = Map::new();
        m.insert("id".into(), Value::from(self.ga4gh_id()));
        m.insert("type".into(), Value::from("CopyNumberChange"));
        m.insert("location".into(), self.location.to_output_value());
        m.insert("copyChange".into(), Value::from(self.copy_change.clone()));
        Value::Object(m)
    }
}

/// A junction between two sequence locations (structural variation). Identifiable —
/// prefix `AJ`. `adjoinedSequences` order is significant (unlike CisPhasedBlock members).
#[derive(Debug, Clone, Deserialize)]
pub struct Adjacency {
    #[serde(rename = "adjoinedSequences")]
    pub adjoined_sequences: Vec<SequenceLocation>,
    /// Optional inserted sequence at the junction.
    #[serde(default)]
    pub linker: Option<State>,
}

impl Adjacency {
    fn digest_value(&self) -> Value {
        let mut m = Map::new();
        // Nested identifiable SequenceLocations are referenced by bare digest, in order.
        let adjoined: Vec<Value> = self
            .adjoined_sequences
            .iter()
            .map(|sl| Value::from(sl.digest()))
            .collect();
        m.insert("adjoinedSequences".into(), Value::Array(adjoined));
        // `linker` is a non-identifiable inline state; always emitted (null when absent).
        m.insert(
            "linker".into(),
            match &self.linker {
                Some(s) => s.digest_value(),
                None => Value::Null,
            },
        );
        m.insert("type".into(), Value::from("Adjacency"));
        Value::Object(m)
    }

    pub fn ga4gh_serialize(&self) -> String {
        serde_json::to_string(&self.digest_value()).expect("serialize Adjacency")
    }

    pub fn digest(&self) -> String {
        sha512t24u(self.ga4gh_serialize().as_bytes())
    }

    pub fn ga4gh_id(&self) -> String {
        format!("ga4gh:AJ.{}", self.digest())
    }

    /// Full VRS 2.0 JSON object with computed identifiers, for NDJSON output.
    pub fn to_output_value(&self) -> Value {
        let adjoined: Vec<Value> = self
            .adjoined_sequences
            .iter()
            .map(|sl| sl.to_output_value())
            .collect();
        let mut m = Map::new();
        m.insert("id".into(), Value::from(self.ga4gh_id()));
        m.insert("type".into(), Value::from("Adjacency"));
        m.insert("adjoinedSequences".into(), Value::Array(adjoined));
        if let Some(s) = &self.linker {
            m.insert("linker".into(), s.digest_value());
        }
        Value::Object(m)
    }
}
