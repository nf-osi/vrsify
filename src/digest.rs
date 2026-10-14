//! GA4GH digest primitive: `sha512t24u`.
//!
//! SHA-512 over the input blob, truncated to the first 24 bytes, then
//! base64url-encoded without padding. Validated against the GA4GH VRS
//! `validation/functions.yaml` golden fixtures (see `tests/`).

use base64::{engine::general_purpose::URL_SAFE_NO_PAD, Engine as _};
use sha2::{Digest, Sha512};

/// Compute the GA4GH `sha512t24u` digest of a byte blob.
pub fn sha512t24u(blob: &[u8]) -> String {
    let mut hasher = Sha512::new();
    hasher.update(blob);
    let full = hasher.finalize();
    URL_SAFE_NO_PAD.encode(&full[..24])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_vectors() {
        // From GA4GH vrs/validation/functions.yaml
        assert_eq!(sha512t24u(b""), "z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc");
        assert_eq!(sha512t24u(b"ACGT"), "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2");
    }
}
