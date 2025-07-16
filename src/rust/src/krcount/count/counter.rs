use std::hash::Hash;

use rustc_hash::FxHashSet as HashSet;

/// Generic trait for counting (either total or unique)
pub(super) trait Countable {
    type Item;
    fn insert(&mut self, item: Self::Item);
    fn count(&self) -> usize;
}

pub(super) struct CountTotal {
    n: usize,
}

impl CountTotal {
    pub(super) fn new() -> Self {
        Self { n: 0 }
    }
}

impl Countable for CountTotal {
    type Item = ();

    fn insert(&mut self, _item: ()) {
        self.n += 1;
    }

    fn count(&self) -> usize {
        self.n
    }
}

pub(super) struct CountUnique<T> {
    set: HashSet<T>,
    n: usize,
}

impl<T> CountUnique<T> {
    #[allow(dead_code)]
    pub(super) fn new() -> Self {
        Self::with_capacity(0)
    }

    pub(super) fn with_capacity(capacity: usize) -> Self {
        Self {
            set: HashSet::with_capacity_and_hasher(capacity, rustc_hash::FxBuildHasher),
            n: 0,
        }
    }
}

impl<T: Eq + Hash> Countable for CountUnique<T> {
    type Item = T;

    fn insert(&mut self, item: T) {
        if self.set.insert(item) {
            self.n += 1;
        }
    }

    fn count(&self) -> usize {
        self.n
    }
}
