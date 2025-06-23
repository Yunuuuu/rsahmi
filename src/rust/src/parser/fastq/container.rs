#[derive(Debug, Default)]
pub struct FastqContainer<'a> {
    id: Option<&'a [u8]>,
    desc: Option<Option<&'a [u8]>>,
    seq: Option<&'a [u8]>,
    sep: Option<&'a [u8]>,
    qual: Option<&'a [u8]>,
}

impl<'a> std::fmt::Display for FastqContainer<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let out = format!(
            "{}{}{}{}",
            match self.id {
                Some(id) => match self.desc {
                    Some(Some(desc)) => format!(
                        "head: {} {}",
                        String::from_utf8_lossy(id),
                        String::from_utf8_lossy(desc)
                    ),
                    _ => format!("head: {}", String::from_utf8_lossy(id)),
                },
                None => String::new(),
            },
            match self.seq {
                Some(seq) =>
                    format!("\nsequence: {}", String::from_utf8_lossy(seq)),
                None => String::new(),
            },
            match self.sep {
                Some(sep) =>
                    format!("\nseparator: {}", String::from_utf8_lossy(sep)),
                None => String::new(),
            },
            match self.qual {
                Some(qual) =>
                    format!("\nquality: {}", String::from_utf8_lossy(qual)),
                None => String::new(),
            }
        );

        f.write_str(&out)
    }
}

impl<'a> FastqContainer<'a> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn reset(&mut self) {
        self.id = None;
        self.desc = None;
        self.seq = None;
        self.sep = None;
        self.qual = None;
    }

    // --- Setters ---
    pub fn set_id(&mut self, id: &'a [u8]) {
        self.id = Some(id);
    }

    pub fn set_desc(&mut self, desc: Option<&'a [u8]>) {
        self.desc = Some(desc);
    }

    pub fn set_seq(&mut self, seq: &'a [u8]) {
        self.seq = Some(seq);
    }

    pub fn set_sep(&mut self, sep: &'a [u8]) {
        self.sep = Some(sep);
    }

    pub fn set_qual(&mut self, qual: &'a [u8]) {
        self.qual = Some(qual);
    }

    // --- Getters ---
    pub fn id(&self) -> Option<&'a [u8]> {
        self.id
    }

    pub fn desc(&self) -> Option<Option<&'a [u8]>> {
        self.desc
    }

    pub fn seq(&self) -> Option<&'a [u8]> {
        self.seq
    }

    pub fn sep(&self) -> Option<&'a [u8]> {
        self.sep
    }

    pub fn qual(&self) -> Option<&'a [u8]> {
        self.qual
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_state_setters_and_getters() {
        let id: &'static [u8] = b"SEQ_ID";
        let desc: &'static [u8] = b"description";
        let seq: &'static [u8] = b"ACGT";
        let sep: &'static [u8] = b"+";
        let qual: &'static [u8] = b"####";

        let mut container = FastqContainer::new();
        container.set_id(id);
        container.set_desc(Some(desc));
        container.set_seq(seq);
        container.set_sep(sep);
        container.set_qual(qual);

        assert_eq!(container.id(), Some(id));
        assert_eq!(container.desc(), Some(Some(desc)));
        assert_eq!(container.seq(), Some(seq));
        assert_eq!(container.sep(), Some(sep));
        assert_eq!(container.qual(), Some(qual));
    }

    #[test]
    fn test_display_with_all_fields() {
        let mut container = FastqContainer::new();
        container.set_id(b"SEQ_ID");
        container.set_desc(Some(b"desc"));
        container.set_seq(b"ACGT");
        container.set_sep(b"+");
        container.set_qual(b"####");

        let formatted = format!("{}", container);
        assert!(formatted.contains("head: SEQ_ID desc"));
        assert!(formatted.contains("sequence: ACGT"));
        assert!(formatted.contains("separator: +"));
        assert!(formatted.contains("quality: ####"));
    }

    #[test]
    fn test_display_with_missing_optional_desc() {
        let mut container = FastqContainer::new();
        container.set_id(b"SEQ_ID");
        container.set_desc(None);
        container.set_seq(b"ACGT");
        container.set_sep(b"+");
        container.set_qual(b"####");

        let formatted = format!("{}", container);
        assert!(formatted.contains("head: SEQ_ID\n"));
        assert!(formatted.contains("sequence: ACGT"));
    }

    #[test]
    fn test_reset() {
        let mut container = FastqContainer::new();
        container.set_id(b"id");
        container.set_desc(Some(b"desc"));
        container.set_seq(b"seq");
        container.set_sep(b"+");
        container.set_qual(b"qual");

        container.reset();
        assert!(container.id().is_none());
        assert!(container.desc().is_none());
        assert!(container.seq().is_none());
        assert!(container.sep().is_none());
        assert!(container.qual().is_none());
    }
}
